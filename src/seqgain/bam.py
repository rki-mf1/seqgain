"""Streaming fragment-depth histograms; overlaps count once per named fragment."""
import csv
from collections import Counter
import pysam
from .io import priority, bam_record_total
from .models import fit_breadth
from .progress import Progress


def read_targets(path, references, lengths):
    known = dict(zip(references, lengths))
    if path:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "contig" not in reader.fieldnames:
                raise ValueError("Reference TSV requires a contig column")
            targets = [{"contig": row["contig"], "group": row.get("group") or "reference"} for row in reader]
        seen = set()
        for t in targets:
            if t["contig"] not in known or t["contig"] in seen:
                raise ValueError(f"Unknown or repeated target contig: {t['contig']}")
            seen.add(t["contig"])
    else:
        targets = [{"contig": c, "group": "reference"} for c in references]
    if not targets:
        raise ValueError("No reference targets selected")
    for t in targets:
        t["length"] = known[t["contig"]]
        if t["length"] <= 0:
            raise ValueError("Reference contigs must have positive lengths")
    return targets


def usable(read, mapq, include_duplicates):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or
                (read.is_duplicate and not include_duplicates) or read.mapping_quality < mapq or
                read.mapping_quality == 255)


def analyze_bam(path, targets_path, fractions, replicates, seed, mapq, baseq,
                include_duplicates=False, max_depth=1000000, progress=None):
    progress = progress or Progress()
    results = []
    display_thresholds = list(range(1, 101))
    with pysam.AlignmentFile(path, "rb") as bam:
        if not bam.has_index():
            raise ValueError("BAM must be coordinate sorted and indexed (samtools index sample.bam)")
        targets = read_targets(targets_path, bam.references, bam.lengths)
        for index, target in enumerate(targets, 1):
            progress.phase(f"Coverage {index}/{len(targets)}: {target['contig']} "
                           f"(MAPQ ≥{mapq}, duplicates {'included' if include_duplicates else 'excluded'})",
                           target["length"], "positions")
            thresholds = display_thresholds
            fit_thresholds = [1, 5, 10, 20]
            hists = [[Counter() for _ in fractions] for _ in range(replicates)]
            visited = 0
            # Pileup streams reference columns, so RAM does not scale with reference length.
            for col in bam.pileup(target["contig"], 0, target["length"], truncate=True,
                                  stepper="nofilter", min_base_quality=0, min_mapping_quality=0,
                                  ignore_overlaps=False, ignore_orphans=False, compute_baq=False,
                                  max_depth=max_depth):
                if col.nsegments >= max_depth:
                    raise ValueError("Pileup depth cap reached; increase --max-pileup-depth to avoid truncation")
                names = set()
                for p in col.pileups:
                    read = p.alignment
                    if p.is_del or p.is_refskip or not usable(read, mapq, include_duplicates):
                        continue
                    q = p.query_position
                    if read.query_qualities is None or read.query_qualities[q] < baseq:
                        continue
                    # RG distinguishes unrelated templates with recycled QNAMEs.
                    rg = str(read.get_tag("RG")) if read.has_tag("RG") else ""
                    names.add((rg, read.query_name or ""))
                visited += 1
                for r in range(replicates):
                    # Hash QNAME so mates (including different contigs) share selection.
                    hashes = sorted(priority(name, seed, r) for name in names)
                    import bisect
                    for j, f in enumerate(fractions):
                        hists[r][j][bisect.bisect_left(hashes, f)] += 1
                progress.update(col.reference_pos + 1)
            progress.update(target["length"])
            progress.phase(f"Fitting coverage forecast {index}/{len(targets)}: {target['contig']}")
            rows = []
            histograms = {}
            for r in range(replicates):
                for j, f in enumerate(fractions):
                    hist = hists[r][j]
                    hist[0] += target["length"] - visited
                    histograms[(r, f)] = dict(hist)
                    rows.append({"replicate": r, "fraction": f,
                                 "mean_fragment_depth": sum(d*n for d,n in hist.items()) / target["length"],
                                 "breadth": {str(d): sum(n for depth,n in hist.items() if depth >= d) /
                                             target["length"] for d in thresholds}})
            full_hist = hists[0][-1]
            result = {**target, "thresholds": thresholds, "curves": rows,
                      "full_depth_histogram": dict(sorted(full_hist.items())),
                      "model": fit_breadth(rows, fit_thresholds, full_hist, histograms, thresholds)}
            results.append(result)
    progress.phase("Combining reference coverage")
    total_length = sum(t["length"] for t in results)
    summary_rows = []
    for replicate in range(replicates):
        for fraction in fractions:
            selected_rows = [next(row for row in target["curves"]
                                  if row["replicate"] == replicate and row["fraction"] == fraction)
                             for target in results]
            summary_rows.append({
                "replicate": replicate,
                "fraction": fraction,
                "mean_fragment_depth": sum(row["mean_fragment_depth"] * target["length"]
                                           for row, target in zip(selected_rows, results)) / total_length,
                "breadth": {str(depth): sum(row["breadth"][str(depth)] * target["length"]
                                            for row, target in zip(selected_rows, results)) / total_length
                            for depth in display_thresholds},
            })
    overall_histogram = Counter()
    for target in results:
        overall_histogram.update({int(depth): count
                                  for depth, count in target["full_depth_histogram"].items()})
    summary = {"label": "All references", "length": total_length,
               "thresholds": display_thresholds, "curves": summary_rows,
               "full_depth_histogram": dict(sorted(overall_histogram.items()))}
    if all(target["model"]["forecast_efforts"] for target in results):
        efforts = results[0]["model"]["forecast_efforts"]
        summary["forecast_efforts"] = efforts
        summary["forecast_breadth"] = {
            str(depth): [sum(target["model"]["forecast_breadth"][str(depth)][i] * target["length"]
                             for target in results) / total_length
                         for i in range(len(efforts))]
            for depth in display_thresholds
        }
    else:
        summary["forecast_efforts"] = []
        summary["forecast_breadth"] = {}
    # Read-level counts explicitly avoid representing these as fragment abundance.
    stats = Counter()
    selected = {r["contig"] for r in results}
    with pysam.AlignmentFile(path, "rb") as bam:
        progress.phase("Counting BAM filter diagnostics", bam_record_total(bam), "records")
        for read in bam.fetch(until_eof=True):
            progress.advance()
            stats["alignment_records"] += 1
            if read.is_secondary:
                stats["filtered_secondary_records"] += 1
                continue
            if read.is_supplementary:
                stats["filtered_supplementary_records"] += 1
                continue
            stats["primary_records"] += 1
            if read.is_unmapped:
                stats["filtered_unmapped_records"] += 1
            elif read.is_qcfail:
                stats["filtered_qcfail_records"] += 1
            elif read.is_duplicate and not include_duplicates:
                stats["filtered_duplicate_records"] += 1
            elif read.mapping_quality == 255:
                stats["filtered_unknown_mapq_records"] += 1
            elif read.mapping_quality < mapq:
                stats["filtered_low_mapq_records"] += 1
            if usable(read, mapq, include_duplicates):
                stats["usable_mapped_records"] += 1
                if read.reference_name in selected:
                    stats["usable_target_records"] += 1
    warnings = ["Coverage is fragment depth: overlapping mates with the same read-group/QNAME count once per position, so it can be lower than read depth.",
                f"Coverage excludes bases below Q{baseq}, records below MAPQ {mapq}, MAPQ 255, QC failures, secondary/supplementary alignments" +
                (", but includes marked duplicates." if include_duplicates else ", and marked duplicates."),
                "Forecast assumes unchanged library, target abundance, mapping and coverage bias.",
                "Absent segments and systematic gaps cannot be distinguished from undersampling without external evidence."]
    if stats["primary_records"] and stats["usable_target_records"] / stats["primary_records"] < .25:
        warnings.append("Fewer than 25% of primary BAM records are usable on the selected references; inspect the filtering counts and target selection when comparing coverage with another tool.")
    return {"targets": results, "summary": summary, "record_counts": dict(stats),
            "coverage_definition": {"unit": "fragment", "overlapping_mates": "count_once",
                                    "denominator": "full BAM-header length of selected references",
                                    "minimum_mapq": mapq, "minimum_baseq": baseq,
                                    "include_marked_duplicates": include_duplicates},
            "warnings": warnings}
