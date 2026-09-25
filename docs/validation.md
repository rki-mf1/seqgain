# Validation

The automated tests cover FASTQ and BAM input handling, read matching, k-mer counting and error diagnostics, deterministic subsampling, fragment-depth coverage and filtering, reference selection, forecast anchoring and consistency, progress output, and generation of an offline HTML report. Synthetic gamma–Poisson data check recovery of known model parameters and retrospective predictions. Report-script interactions are checked with Node.js; the Conda test environment includes it.

GitHub Actions builds and installs the Conda package in a Python 3.12 environment, runs the tests and CLI smoke checks, and saves a generated demonstration report as an expiring artifact. It also checks for merge-conflict markers.

These synthetic checks verify implementation behavior, not the predictive accuracy of forecasts on biological samples. Broad real-data calibration and production-scale benchmarks have not been documented. Treat forecasts as exploratory, particularly beyond the tested 2× retrospective horizon; compare them with withheld or independent sequencing data before using them to guide decisions.
