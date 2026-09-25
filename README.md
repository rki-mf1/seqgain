# SeqGain

## Abstract

SeqGain helps explore what additional sequencing might reveal. It measures how many new short sequence patterns (k-mers) appear as more reads are examined, and what fraction of each reference reaches a chosen coverage depth. An interactive, offline HTML report shows the measurements alongside projections for more sequencing.

Use a BAM file for both analyses, or FASTQ files for k-mer discovery; reference coverage requires a BAM file. SeqGain is particularly useful for segmented genomes, plasmids and fragmented assemblies, where an overall average could hide a poorly covered reference. Projections assume that future reads resemble the current data. They are exploratory estimates, not validated stopping rules or sequencing recommendations.

Project concept by **Kathrin Trappe** and **Simon H. Tausch**.

SeqGain is research software, not clinically validated. It does not establish assembly completeness, resolve repeats or infer plasmid copy number. The report runs locally without uploads, accounts or telemetry.

## Installation

Build and install from a checkout of this repository:

```bash
conda env create -f environment.yml
conda activate seqgain
conda-build conda-recipe --output-folder ./conda-channel --override-channels -c conda-forge -c bioconda
conda install --strict-channel-priority --override-channels -c "file://$PWD/conda-channel" -c conda-forge -c bioconda seqgain=0.3.0
seqgain --help
```

Core dependencies are NumPy, SciPy, pysam/HTSlib and Plotly. The source-development environment also includes Minimap2 and samtools for optional alignment preparation; they are not required by the SeqGain package itself.

## Try the included demo

Generate the deterministic synthetic inputs and a report in new directories:

```bash
seqgain-demo --outdir demo-input
seqgain \
  --reads1 demo-input/reads.fastq \
  --bam demo-input/sample.bam \
  --targets demo-input/targets.tsv \
  --sample 'Segmented virus demo' \
  --outdir demo-report
```

Open `demo-report/report.html` in your browser. The synthetic demo has two segments at different depths; segment 2 is deliberately less covered at a 10× minimum depth.

## Use your data

Inputs: four-line Phred+33 FASTQ (plain text or `.gz`) and/or coordinate-sorted, indexed BAM. CRAM, FASTA reads, interleaved FASTQ and remote URLs are unsupported. Output directories must be new or empty.

Paired-end reads with a matching reference alignment:

```bash
seqgain --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \
  --bam sample.sorted.bam \
  --sample isolate_A --outdir results/isolate_A
```

Reference-free only:

```bash
seqgain --reads1 sample.fastq.gz --k 31 --outdir results/kmer
```

Both coverage and k-mer discovery from BAM alone (default):

```bash
seqgain --bam sample.sorted.bam --outdir results/reference
```

K-mer input defaults to BAM when present, otherwise FASTQ. Use `--kmer-source fastq` to prefer supplied FASTQs or `--kmer-source none` for coverage only.

BAM k-mer discovery uses all stored primary reads, including unmapped, duplicate, QC-failed and low-MAPQ reads; only base-quality filtering applies. If the BAM lacks full sequences or qualities (for example, because reads were hard-clipped), use FASTQ instead. Reads removed before BAM creation cannot be recovered.

For joint interpretation, BAM and FASTQ should represent the same sequencing effort. When both are supplied, SeqGain compares read counts, names, sequences and qualities; a mismatch produces a warning in the report. This check does not prove biological sample identity.

### Explore alignment filters in HTML

Add `--explore-filters` to a BAM run to enable minimum-MAPQ and marked-duplicate controls in the report:

```bash
seqgain --bam sample.sorted.bam --explore-filters --outdir results/filter-exploration
```

This precomputes MAPQ **0, 20, 30** (plus a custom `--mapq` value) with duplicates included and excluded. Switching filters updates coverage results instantly in the offline report; k-mer results do not change. MAPQ 255 is always excluded.

This requires **6 coverage analyses (8 with a custom MAPQ)**, so runtime and report size increase. Without the flag, the report uses the fixed CLI filters. See [Outputs](#outputs) for exported variants.

### Forecast range

The report initially shows 0–2× the current sequencing effort. One slider extends the plot range to 100×; another selects a scenario. Values below 1× come from observed subsamples, 1× is the full input, and values above 1× are forecasts. The 2× default is only a display choice: longer-range forecasts are available but not validated by the retrospective 2× checks. A separate 1–100× minimum-depth slider changes which reference-breadth threshold is shown, not the model fit.

### Select reference segments

Metagenome with a focal segmented virus:

```bash
seqgain --reads1 meta_R1.fastq.gz --reads2 meta_R2.fastq.gz \
  --bam meta.sorted.bam --targets targets.tsv \
  --sample metagenome_A --outdir results/metagenome_A
```

K-mers describe the **whole supplied read set**, while coverage describes only selected references. For target-specific discovery, supply consistently target-selected FASTQs and BAM, and interpret effort relative to that subset. No taxonomic classifier is run.

The targets TSV selects BAM header references; `contig` is required and `group` is an optional label:

```tsv
contig	group
segment_1	virus
segment_2	virus
```

Without this file, all BAM references are included. Each reference has its own forecast; the overall view is length-weighted.

### Analysis options

Run `seqgain --help` for all options and defaults. Key controls are:

| Option | Meaning |
|---|---|
| `--baseq 20` | Minimum Phred base quality. A k-mer is retained only when all its bases pass; a BAM position is covered only by passing bases. |
| `--mapq 30` | Minimum BAM alignment mapping quality. Lowering it admits more ambiguously mapped alignments. MAPQ 255 means “mapping quality unavailable” and is always excluded, not treated as a very confident score. |
| `--include-duplicates` | Include BAM records marked as PCR/optical duplicates. They are excluded by default; unmarked duplicates cannot be detected reliably here. |
| `--explore-filters` | Precompute MAPQ/duplicate variants for interactive coverage filtering; see overhead above. Requires BAM. |
| `--quiet` | Hide the two-second terminal progress updates, but retain warnings and the report path. |
| `--targets FILE` | Analyze only listed references and optionally label their groups. Without it, every BAM header reference is included. |
| `--fractions ...` | Nested sampling fractions for discovery and coverage curves. Default: `0.05,0.1,0.2,0.25,0.4,0.5,0.6,0.8,1`. Must contain at least four distinct values in `(0,1]`, including `1` (the full input). |
| `--replicates 3` / `--seed 42` | Number and seed of deterministic hash-subsampling replicates. These show thinning variability, not biological replication. |
| `--k 31` | Canonical k-mer length for BAM or FASTQ discovery. |
| `--max-pileup-depth 1000000` | Per-position safety cap. SeqGain fails if it is reached instead of silently truncating deep coverage. |

SeqGain reports **fragment coverage**, not conventional read coverage: overlapping mates with the same read-group and QNAME count once at a position. Its denominator is the full selected reference length, including zero-coverage bases. When comparing with `samtools depth` or another manual calculation, use matching MAPQ/base-quality, duplicate, overlap and zero-depth settings. The report records filtering counts and the exact coverage definition under provenance to help diagnose discrepancies.

### Prepare an alignment if needed

Example for Illumina short paired reads (a Bash pipeline):

```bash
set -o pipefail
minimap2 -t 8 -ax sr reference.fa sample_R1.fastq.gz sample_R2.fastq.gz \
  | samtools sort -@ 4 -o sample.sorted.bam
samtools index sample.sorted.bam
```

Choose a mapper/preset appropriate for your technology. In metagenomes, competitive mapping against background/decoy sequences can help assess ambiguous matches. Host depletion, adapter trimming, duplicate marking, alignment QC and reference selection remain upstream responsibilities; unmarked PCR duplicates cannot be identified from QNAME alone.

## Outputs

| File | Contents |
|---|---|
| `report.html` | Standalone interactive report with bundled JavaScript |
| `results.json` | Metrics, histograms, forecasts, diagnostics, parameters and provenance; default coverage under `bam`, optional alternatives under `bam_filter_variants` |
| `curves.tsv` | Observed replicate curves, including `all_references`; only the CLI-selected coverage filters |
| `kmers.sqlite` | Optional exact count database with `--keep-kmer-db` |

The HTML JSON download includes all filter variants, not just the displayed setting. Provenance records versions, command and input paths/sizes/mtimes—not file-content checksums. The separate BAM/FASTQ comparison uses read-content fingerprints as described above. Mapping counts represent BAM records, not template abundance; a mapped-only BAM cannot establish the original mapping fraction. Failures leave partial files for diagnosis; rerun into a new directory rather than overwriting results.

## Interpretation and limits

Compare reference breadth and k-mer novelty at the same effort. Marginal k-mer gain is also shown for the next 0.1× baseline effort. SeqGain does not make an automated sequencing recommendation.

K-mer richness is neither molecular library complexity nor metagenome completeness. Errors, contamination, uneven abundances, repeats, k and quality thresholds affect it. High singleton fractions trigger a warning. Recurrent k-mers are reported separately; singletons remain in the discovery curve.

The report also calculates a Phred-based upper estimate of how many k-mers could arise from base-call errors. This diagnostic does not change the counts or Chao2 forecast; it cannot detect every source of error or distinguish errors from genuinely rare sequence. See [Methods](docs/methods.md) for its assumptions and calculation.

Reference forecasts start at the measured breadth and use the observed depth distribution to project additional coverage. Breadth at a higher depth threshold cannot exceed breadth at a lower one. Poor fits and retrospective prediction errors are flagged. Forecasts assume unchanged sample composition, mapping and coverage bias; they are **not guarantees or confidence intervals**. The three default subsampling replicates show variation when reads are thinned, not uncertainty about unseen sequence. References without usable coverage have no forecast.

### Related work

SeqGain was developed with awareness of [preseq](https://preseq.readthedocs.io/en/latest/) and [preseqR](https://github.com/smithlabcode/preseqR). We acknowledge the contributions of their authors to library-complexity, genome-coverage and higher-depth coverage prediction, including [Daley & Smith (2013)](https://doi.org/10.1038/nmeth.2375), [Daley & Smith (2014)](https://doi.org/10.1093/bioinformatics/btu540) and [Deng et al. (2020)](https://doi.org/10.1089/cmb.2019.0264). Predicting the benefit of additional sequencing is not a new concept introduced by SeqGain.

Our intended contribution is a complementary workflow: per-reference depth/breadth scenarios alongside whole-input k-mer discovery, with interactive alignment-filter comparisons and diagnostics in one offline report. SeqGain is not a replacement for preseq and makes no claim of superior predictive accuracy.

[Nonpareil](https://nonpareil.readthedocs.io/en/latest/) addresses abundance-weighted metagenomic coverage. It is also a separate tool, not a SeqGain dependency or integrated backend.

### Scale

SQLite bounds k-mer RAM, but disk usage grows with distinct k-mers and replicate count. Exact Python enumeration and database updates can be slow on large metagenomes. Streaming BAM pileups suit viruses, microbial isolates and focal references; per-base Python processing is not optimized for deep whole-human-genome data. Peak pileup memory grows with local depth; the safety cap fails explicitly instead of truncating coverage. Production-scale performance is not established.

See [Methods](docs/methods.md) for equations and assumptions and [Validation](docs/validation.md) for verification details.

### Development note

SeqGain was developed with substantial assistance from OpenAI Codex for code implementation, debugging, test development and documentation. The human authors defined the scientific objectives, directed development, and reviewed the implementation and its assumptions. The authors retain responsibility for the software, its interpretation and its limitations. Automated tests and validation checks are documented separately; AI assistance does not constitute independent scientific validation.

## License

SeqGain is licensed under the GNU General Public License v3.0 only (`GPL-3.0-only`). See [LICENSE](LICENSE) for the full text.
