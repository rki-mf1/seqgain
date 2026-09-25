# Changelog

## 0.3.1 — SeqGain release

* Publish the renamed SeqGain package with GitHub Actions and GitHub source archives.
* Reject FASTQ records without a read name with a clear validation error.
* Reject runs that disable k-mer discovery without providing a BAM for coverage.
* Correct the Bioconda release and submission instructions.

## 0.3.0 — initial Bioconda-targeted release

* Exact canonical k-mer discovery from BAM or FASTQ, with quality diagnostics and incidence-based extrapolation.
* Fragment-depth reference coverage, per-segment and overall summaries, anchored effort scenarios, and retrospective checks.
* Standalone interactive HTML report with optional precomputed alignment-filter comparisons.
* Conda build recipe, automated tests, and documented scientific assumptions and limitations.
