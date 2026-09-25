"""Exact disk-backed canonical k-mer incidence, counted once per fragment."""
import json
import math
import sqlite3
from collections import Counter, deque
from pathlib import Path
from .io import fragments, priority
from .progress import Progress

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical_kmer_evidence(reads, k, min_baseq):
    """Return canonical k-mers and a conservative Phred error probability.

    Values are upper bounds on the probability that every occurrence of the
    word in this fragment contains a base-calling error. For repeated words we
    retain the best-quality occurrence; this remains conservative when windows
    overlap and their base-call errors are dependent.
    """
    evidence = {}
    for _, sequence, quality in reads:
        window = deque()
        log_correct = 0.0
        certain_error_bases = 0
        for i, (base, q) in enumerate(zip(sequence, quality)):
            score = ord(q) - 33
            if base not in "ACGT" or score < min_baseq:
                window.clear()
                log_correct = 0.0
                certain_error_bases = 0
                continue
            error_probability = 10 ** (-score / 10)
            value = None if error_probability >= 1 else math.log1p(-error_probability)
            window.append(value)
            if value is None:
                certain_error_bases += 1
            else:
                log_correct += value
            if len(window) > k:
                old = window.popleft()
                if old is None:
                    certain_error_bases -= 1
                else:
                    log_correct -= old
            if len(window) == k:
                word = sequence[i-k+1:i+1]
                word = min(word, word.translate(COMPLEMENT)[::-1])
                p_wrong = 1.0 if certain_error_bases else -math.expm1(log_correct)
                evidence[word] = min(evidence.get(word, 1.0), p_wrong)
    return evidence


def canonical_kmers(reads, k, min_baseq):
    return set(canonical_kmer_evidence(reads, k, min_baseq))


def discovery_forecast(n, observed, q1, q2, multipliers):
    # Bias-corrected incidence Chao2 lower-bound unseen richness.
    unseen = ((n - 1) / n) * q1 * (q1 - 1) / (2 * (q2 + 1)) if n else 0
    rate = q1 / (n * unseen + q1) if unseen and q1 else 0
    return {"unseen_lower_bound": unseen, "asymptote_lower_bound": observed + unseen,
            "points": [{"effort": float(t), "distinct": observed + unseen *
                        (-math.expm1(n * (t - 1) * math.log1p(-rate))) if 0 < rate < 1 else observed}
                       for t in multipliers]}


def analyze_kmers(r1, r2, outdir, fractions, replicates, seed, k=31, min_baseq=20, fragment_stream=None,
                  progress=None):
    progress = progress or Progress()
    progress.phase("Counting k-mer reads", unit="reads")
    dbpath = Path(outdir) / "kmers.sqlite"
    if dbpath.exists():
        raise ValueError(f"Refusing to overwrite {dbpath}")
    conn = sqlite3.connect(dbpath)
    cols = [[f"h{i}_{rank}" for rank in range(1, 4)] for i in range(replicates)]
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-32768")
    conn.execute("CREATE TABLE names (name TEXT PRIMARY KEY) WITHOUT ROWID")
    rank_fields = [c for group in cols for c in group]
    conn.execute("CREATE TABLE kmers (k TEXT PRIMARY KEY, n INTEGER NOT NULL, "
                 "false_probability REAL NOT NULL," +
                 ",".join(f"{c} REAL NOT NULL" for c in rank_fields) + ") WITHOUT ROWID")
    updates = ["n=n+1", "false_probability=false_probability*excluded.false_probability"]
    for first, second, third in cols:
        updates.extend([
            f"{third}=CASE WHEN excluded.{first}<{first} THEN {second} "
            f"WHEN excluded.{first}<{second} THEN {second} ELSE min({third},excluded.{first}) END",
            f"{second}=CASE WHEN excluded.{first}<{first} THEN {first} ELSE min({second},excluded.{first}) END",
            f"{first}=min({first},excluded.{first})",
        ])
    sql = "INSERT INTO kmers VALUES (" + ",".join("?" for _ in range(3+3*replicates)) + ") " + \
          "ON CONFLICT(k) DO UPDATE SET " + ",".join(updates)
    total = bases = valid_fragments = 0
    expected_error_incidences = 0.0
    sampled = [[0]*len(fractions) for _ in range(replicates)]
    try:
        for name, reads in (fragment_stream if fragment_stream is not None else fragments(r1, r2)):
            try:
                identity = json.dumps(name, separators=(",", ":")) if isinstance(name, tuple) else name
                conn.execute("INSERT INTO names VALUES (?)", (identity,))
            except sqlite3.IntegrityError as exc:
                raise ValueError(f"Duplicate fragment name in read input: {name}") from exc
            total += 1
            bases += sum(len(read[1]) for read in reads)
            hashes = [priority(name, seed, r) for r in range(replicates)]
            for r, h in enumerate(hashes):
                for j, f in enumerate(fractions):
                    sampled[r][j] += h < f
            evidence = canonical_kmer_evidence(reads, k, min_baseq)
            valid_fragments += bool(evidence)
            expected_error_incidences += sum(evidence.values())
            rows = []
            for word, p_wrong in evidence.items():
                ranks = [value for h in hashes for value in (h, 2.0, 2.0)]
                rows.append((word, 1, p_wrong, *ranks))
            conn.executemany(sql, rows)
            progress.advance(len(reads))
            if total % 1000 == 0:
                conn.commit()
        conn.commit()
        if not total:
            raise ValueError("Read input contains no fragments")
        progress.phase("Summarizing k-mer spectrum")
        spectrum = dict(conn.execute("SELECT n,count(*) FROM kmers GROUP BY n"))
        distinct = sum(spectrum.values())
        if not distinct:
            raise ValueError("No valid k-mers; check k, read length and quality threshold")
        curves = []
        progress.phase("Calculating k-mer subsamples", replicates * len(fractions), "subsamples")
        for r, (first, second, third) in enumerate(cols):
            for j, f in enumerate(fractions):
                count = conn.execute(f"SELECT count(*) FROM kmers WHERE {first} < ?", (f,)).fetchone()[0]
                q1_sub = conn.execute(
                    f"SELECT count(*) FROM kmers WHERE {first} < ? AND {second} >= ?", (f, f)).fetchone()[0]
                q2_sub = conn.execute(
                    f"SELECT count(*) FROM kmers WHERE {second} < ? AND {third} >= ?", (f, f)).fetchone()[0]
                curves.append({"replicate": r, "fraction": f, "fragments": sampled[r][j],
                               "distinct": count, "singletons": q1_sub, "doubletons": q2_sub})
                progress.advance()
        progress.phase("Fitting k-mer forecast and diagnostics")
        q1, q2 = spectrum.get(1, 0), spectrum.get(2, 0)
        from .models import effort_grid
        forecast = discovery_forecast(total, distinct, q1, q2, effort_grid())
        expected_false_distinct = conn.execute("SELECT coalesce(sum(false_probability),0) FROM kmers").fetchone()[0]
        expected_false_singletons = conn.execute(
            "SELECT coalesce(sum(false_probability),0) FROM kmers WHERE n=1").fetchone()[0]
        curve_lookup = {(row["replicate"], round(row["fraction"], 12)): row for row in curves}
        backtests = []
        for row in curves:
            target_fraction = round(2 * row["fraction"], 12)
            target = curve_lookup.get((row["replicate"], target_fraction))
            training_count = sum(f <= row["fraction"] for f in fractions)
            if target is None or training_count < 4 or not row["fragments"]:
                continue
            estimate = discovery_forecast(row["fragments"], row["distinct"], row["singletons"],
                                          row["doubletons"], [2])["points"][0]["distinct"]
            backtests.append({"replicate": row["replicate"], "training_fraction": row["fraction"],
                              "validation_fraction": target["fraction"], "predicted_distinct": estimate,
                              "observed_distinct": target["distinct"],
                              "relative_error": (estimate-target["distinct"]) / target["distinct"]})
        warnings = ["K-mer richness is not genome breadth or molecular library complexity.",
                    "Overlapping k-mers are dependent; forecasts are exploratory incidence-model estimates.",
                    "Phred-derived error values assume calibrated base qualities and independent errors across fragments."]
        if q1 / distinct > .5:
            warnings.append("Over half of distinct k-mers are singletons: errors or rare taxa may dominate extrapolation.")
        if q2 == 0:
            warnings.append("No doubletons: unseen-richness estimate is poorly constrained.")
        if backtests and max(abs(test["relative_error"]) for test in backtests) > .1:
            warnings.append("A retrospective 2× k-mer forecast missed observed richness by more than 10%.")
        return {"k": k, "fragments": total, "bases": bases, "valid_fragments": valid_fragments,
                "distinct": distinct, "recurrent_distinct": distinct-q1, "singletons": q1,
                "doubletons": q2, "incidences": sum(n*c for n,c in spectrum.items()),
                "phred_error": {"expected_erroneous_incidences_upper": expected_error_incidences,
                    "expected_false_distinct_upper": expected_false_distinct,
                    "expected_false_singletons_upper": expected_false_singletons,
                    "singleton_fraction_upper": expected_false_singletons/q1 if q1 else 0,
                    "method": "Per-window Phred probability; best occurrence per fragment; independent fragments"},
                "spectrum": spectrum, "curves": curves, "forecast": forecast,
                "retrospective_backtests": backtests, "warnings": warnings}
    finally:
        conn.close()
