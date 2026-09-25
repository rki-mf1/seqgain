import argparse
import csv
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from . import __version__
from .io import fingerprint


def parser():
    p = argparse.ArgumentParser(
        description="SeqGain: analyze FASTQ, indexed BAM, or both",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=__version__)
    p.add_argument("--reads1", type=Path, help="Single-end FASTQ, or mate 1 of paired-end FASTQ (.gz accepted)")
    p.add_argument("--reads2", type=Path, help="Mate 2 FASTQ; requires --reads1 with matching read names")
    p.add_argument("--bam", type=Path, help="Coordinate-sorted, indexed BAM used for reference coverage")
    p.add_argument("--kmer-source", choices=["auto", "bam", "fastq", "none"], default="auto",
                   help="K-mer input: auto prefers BAM, otherwise FASTQ; none disables discovery")
    p.add_argument("--targets", type=Path, help="Optional TSV with contig and optional group columns; otherwise use all BAM references")
    p.add_argument("--outdir", type=Path, required=True, help="New or empty output directory")
    p.add_argument("--sample", default="Sample", help="Sample label displayed in the report")
    advanced = p.add_argument_group("Advanced sampling and storage")
    advanced.add_argument("--fractions", default="0.05,0.1,0.2,0.25,0.4,0.5,0.6,0.8,1",
                   help="Comma-separated nested subsampling fractions, including 1")
    advanced.add_argument("--replicates", type=int, default=3, help="Deterministic hash-subsampling replicates (1-10)")
    advanced.add_argument("--seed", type=int, default=42, help="Seed for deterministic fragment subsampling")
    p.add_argument("--k", type=int, default=31, help="Canonical k-mer length (3-127)")
    p.add_argument("--baseq", type=int, default=20,
                   help="Minimum Phred base quality for FASTQ k-mers and BAM-covered bases (0-93)")
    p.add_argument("--mapq", type=int, default=30,
                   help="Minimum BAM alignment mapping quality (0-254); lower values admit more ambiguous mappings, while MAPQ 255 is always excluded as unknown")
    p.add_argument("--include-duplicates", action="store_true",
                   help="Include BAM records marked as PCR/optical duplicates")
    p.add_argument("--explore-filters", action="store_true",
                   help="Precompute HTML coverage controls for MAPQ 0, 20, 30 and --mapq, each with/without marked duplicates; runs coverage 6–8 times and enlarges the report (k-mers counted once)")
    p.add_argument("--quiet", action="store_true", help="Disable terminal progress; warnings and final report path remain visible")
    advanced.add_argument("--max-pileup-depth", type=int, default=1000000,
                   help="Safety cap per reference position; reaching it fails instead of truncating silently")
    advanced.add_argument("--keep-kmer-db", action="store_true", help="Retain the exact intermediate SQLite k-mer database")
    return p


def demo_main(argv=None):
    p = argparse.ArgumentParser(prog="seqgain-demo", description="Create deterministic synthetic SeqGain inputs",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--outdir", type=Path, required=True, help="New or empty directory for the demo FASTQ, BAM and targets")
    a = p.parse_args(argv)
    try:
        from .demo import make_demo
        make_demo(a.outdir)
    except (ValueError, OSError) as exc:
        p.exit(2, f"seqgain-demo: {exc}\n")


def validate(a):
    if not a.reads1 and not a.bam:
        raise ValueError("Provide --reads1, --bam, or both")
    if a.reads2 and not a.reads1:
        raise ValueError("--reads2 requires --reads1")
    if a.targets and not a.bam:
        raise ValueError("--targets requires --bam")
    if a.explore_filters and not a.bam:
        raise ValueError("--explore-filters requires --bam")
    if a.kmer_source == "bam" and not a.bam or a.kmer_source == "fastq" and not a.reads1:
        raise ValueError("Selected --kmer-source requires the corresponding input")
    a.fractions = sorted(set(float(x) for x in a.fractions.split(",")))
    if len(a.fractions) < 4 or a.fractions[-1] != 1 or not all(0 < x <= 1 for x in a.fractions):
        raise ValueError("Use ≥4 distinct fractions in (0,1], including 1")
    if not 1 <= a.replicates <= 10 or not 3 <= a.k <= 127:
        raise ValueError("Replicates must be 1–10; k must be 3–127")
    if not 0 <= a.baseq <= 93 or not 0 <= a.mapq <= 254:
        raise ValueError("Phred+33 base quality must be 0–93; MAPQ must be 0–254")
    if a.max_pileup_depth < 1:
        raise ValueError("Pileup cap must be positive")
    for path in [a.reads1, a.reads2, a.bam, a.targets]:
        if path and not path.is_file():
            raise ValueError(f"Input not found: {path}")
    if a.outdir.exists() and any(a.outdir.iterdir()):
        raise ValueError("Output directory must be absent or empty; existing results are not overwritten")


def analyze(a):
    from .progress import Progress
    with Progress(enabled=not a.quiet) as progress:
        _analyze(a, progress)


def _analyze(a, progress):
    from .kmers import analyze_kmers
    from .bam import analyze_bam
    from .report import write_report
    validate(a)
    a.outdir.mkdir(parents=True, exist_ok=True)
    parameters = {k: str(v) if isinstance(v, Path) else v for k,v in vars(a).items()}
    result = {"schema_version": 4, "sample": a.sample, "parameters": parameters,
              "provenance": {"seqgain": __version__, "python": platform.python_version(),
                  "utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
                  "dependencies": {n: importlib.metadata.version(n) for n in ["numpy", "scipy", "pysam", "plotly"]},
                  "inputs": [fingerprint(p) for p in [a.reads1,a.reads2,a.bam,a.targets] if p],
                  "input_identity": "Path, size and mtime; not a cryptographic content checksum"},
              "kmers": None, "bam": None, "bam_filter_variants": []}
    source = a.kmer_source if a.kmer_source != "auto" else "bam" if a.bam else "fastq"
    if source != "none":
        from .io import bam_fragments
        stream = bam_fragments(a.bam, a.outdir, progress) if source == "bam" else None
        result["kmers"] = analyze_kmers(a.reads1, a.reads2, a.outdir, a.fractions, a.replicates, a.seed, a.k, a.baseq, stream, progress)
        result["kmers"]["source"] = source
        if source == "bam":
            result["kmers"]["warnings"].append("BAM k-mers use all stored primary reads, including unmapped reads, duplicates and low-MAPQ reads; base quality filtering applies. Reads removed upstream cannot be recovered; mapped-only BAM describes only the retained subset.")
    if a.bam:
        result["bam"] = analyze_bam(a.bam, a.targets, a.fractions, a.replicates, a.seed, a.mapq, a.baseq,
                                   a.include_duplicates, a.max_pileup_depth, progress)
        if a.explore_filters:
            for minimum_mapq in sorted({0, 20, 30, a.mapq}):
                for duplicates in (False, True):
                    if (minimum_mapq, duplicates) == (a.mapq, a.include_duplicates):
                        continue
                    result["bam_filter_variants"].append(analyze_bam(
                        a.bam, a.targets, a.fractions, a.replicates, a.seed, minimum_mapq, a.baseq,
                        duplicates, a.max_pileup_depth, progress))
    if result["kmers"]:
        result["kmers"]["warnings"].append("K-mers describe all supplied reads, not just selected reference contigs.")
    if a.reads1 and a.bam:
        from .io import compare_bam_fastq
        comparison = compare_bam_fastq(a.bam, a.reads1, a.reads2, progress)
        result["provenance"]["bam_fastq_comparison"] = comparison
        if not comparison["match"]:
            warning = "BAM/FASTQ mismatch: primary read counts or name/sequence/quality fingerprints differ. Joint effort comparisons may not describe the same library."
            for coverage in [result["bam"], *result["bam_filter_variants"]]:
                coverage["warnings"].append(warning)
            progress.finish()
            print("WARNING: " + warning, file=sys.stderr)
    progress.phase("Writing JSON, tables and offline HTML")
    (a.outdir / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    with (a.outdir / "curves.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["method", "target", "replicate", "fraction", "metric", "value"])
        if result["kmers"]:
            for r in result["kmers"]["curves"]:
                writer.writerow(["kmer", "all_reads", r["replicate"], r["fraction"], "distinct", r["distinct"]])
        if result["bam"]:
            for t in [result["bam"]["summary"], *result["bam"]["targets"]]:
                for r in t["curves"]:
                    target = t.get("contig", "all_references")
                    writer.writerow(["reference", target, r["replicate"], r["fraction"],
                                     "mean_fragment_depth", r["mean_fragment_depth"]])
                    for depth,value in r["breadth"].items():
                        writer.writerow(["reference", target, r["replicate"], r["fraction"],
                                         "breadth_ge_"+depth, value])
    write_report(result, a.outdir / "report.html")
    if result["kmers"] and not a.keep_kmer_db:
        (a.outdir / "kmers.sqlite").unlink()
    progress.finish()
    print(a.outdir / "report.html")


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    # Compatibility with 0.2 command lines; analysis is now the default and the
    # demo has its own seqgain-demo entry point.
    if args[:1] == ["analyze"]:
        args.pop(0)
    elif args[:1] == ["demo"]:
        return demo_main(args[1:])
    p = parser()
    a = p.parse_args(args)
    try:
        analyze(a)
    except (ValueError, OSError) as exc:
        p.exit(2, f"seqgain: {exc}\n")
