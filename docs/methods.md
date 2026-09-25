# Methods and scientific scope

## Sampling unit and deterministic thinning

FASTQ mates are validated and counted as one fragment. Duplicate normalized fragment identifiers in FASTQ are rejected. Canonical k-mers passing every base's Phred+33 threshold are counted once per fragment, even if repeated within a read or shared by the two mates. Reverse complements are collapsed; ambiguity breaks a k-mer window.

A BLAKE2b hash of seed, replicate and normalized QNAME gives a uniform priority in [0,1). A fragment belongs to a subsample when its priority is below the requested fraction. Samples are nested within a replicate; sample counts are binomial, not fixed. Mates share priority. Different replicates share input data and must not be treated as independent biological replicates. Fractions include the full sample. FASTQ names with `/1` or `/2` are normalized. Unique QNAMEs across the original library are expected; reused names across BAM read groups share selection but are counted as separate fragments at a reference position.

The SQLite table stores total incidence and the smallest fragment priority per replicate for every k-mer. A k-mer is observed at fraction f exactly when its minimum priority is below f. This gives exact discovery counts without retaining complete read sets.

The three smallest priorities are retained for each replicate. They are sufficient to recover singleton and doubleton counts at every configured fraction, allowing a smaller subsample to predict its matching 2× subsample retrospectively. Default fractions include 0.25→0.5 and 0.5→1 comparisons. These checks reuse nested reads and therefore measure forecast calibration on this sample, not independent biological replication.

## Phred-derived k-mer error diagnostic

For every retained k-mer occurrence with base qualities `Q_i`, convert each score to an estimated error probability `e_i=10^(-Q_i/10)`. Assuming calibrated qualities, the probability that all calls in that occurrence are correct is

```
P(correct occurrence) = product_i (1 - e_i),
P(wrong occurrence) = 1 - P(correct occurrence).
```

If the same canonical word occurs more than once in a fragment, SeqGain retains the lowest error probability. This is a conservative upper bound on the probability that all occurrences in that fragment are wrong, including when overlapping windows share base calls. Across different fragments, these probabilities are multiplied under an independence assumption. Summing the resulting values gives upper estimates for false distinct and false singleton k-mers. Summing the per-fragment values gives the expected erroneous incidence burden.

The diagnostic does not remove k-mers or modify Chao2 inputs. Base qualities may be miscalibrated, and PCR errors, contamination, index misassignment and systematic errors may receive high quality. Conversely, a genuine sequence from a rare organism can be a singleton. The estimate is therefore an error baseline rather than a classifier.

## K-mer extrapolation

Let N be the number of input fragments, S the number of distinct observed k-mers, and Q1 and Q2 the numbers observed in exactly one and two fragments. The bias-corrected incidence Chao2 estimate of unseen richness is

```
Q0 = (N - 1)/N × Q1(Q1 - 1) / [2(Q2 + 1)].
```

For m=N(t−1) additional fragments, the implemented incidence extrapolation is

```
S(t) = S + Q0 × [1 − (1 − Q1/(N Q0 + Q1))^m],  1 ≤ t ≤ 100.
```

When Q0=0 the curve stays at S. S+Q0 is a lower-bound richness estimate under the incidence sampling model, not a known population total. Sparse/heterogeneous data can underestimate unseen richness; sequencing errors can inflate it. Shared/overlapping k-mers violate independence assumptions, so this application reports an exploratory point scenario and does not attach confidence intervals. Per-fragment incidence helps avoid treating a homopolymer or overlapping mates as many independent detections, but it does not eliminate dependence between k-mers. All fragments, including those with no passing windows, count toward N and effort.

This is sequence-feature discovery, not the number of original molecules. For molecular library complexity, use a method with appropriate molecule/UMI or positional duplicate definitions.

## Reference breadth

For every BAM reference column, exclude unmapped, QC-failed, secondary, supplementary, below-MAPQ and unknown-MAPQ (255) alignments; marked duplicates are excluded unless requested. Exclude low-quality or missing-quality bases. CIGAR deletions and reference skips do not cover bases. Soft clips and insertions do not cover extra reference bases. No BAQ modification is applied. Overlapping mates count once per `(RG,QNAME)` at each position, even if their sequences disagree; this is fragment depth, not a consensus correctness metric. Orphaned primary alignments can contribute.

Reference length includes all header positions, including zero-depth positions. At each fraction and replicate the histogram H(d) defines

```
mean fragment depth = sum_d d H(d) / reference_length
breadth at threshold q = sum_{d >= q} H(d) / reference_length.
```

Each segment has independent model parameters; there is no equal-abundance assumption or inferred ploidy.

## Gamma–Poisson reference forecast

Accessible fraction a has site-specific intensities drawn from a gamma distribution with mean μ and shape α. Remaining fraction 1−a has zero intensity. At relative sequencing effort t, depth on accessible sites is negative binomial with size α and success probability α/(α+μt). Thus

```
B_q(t) = a × Pr[NB(size=α, p=α/(α+μt)) >= q].
```

Breadth is recorded for every integer depth from 1× through 100×. The three model parameters are fitted by bounded nonlinear least squares to mean replicate breadth across fractions and thresholds 1,5,10,20, independent of the display slider. Bounds are a∈[10⁻⁶,1], μ∈[10⁻⁶,10⁶], α∈[0.03,10⁴]; three dispersion starts reduce local-fit sensitivity. All observations receive equal weight. Correlated fractions/thresholds mean residuals are diagnostics, not inferential independent errors.

A separate fit excluding fraction 1 predicts the full-effort breadth at each threshold. Report its maximum absolute error and the all-data RMSE. Flag RMSE>0.03, holdout error>0.05, optimizer failure, Jacobian condition number>10⁶, or certain boundary estimates.

Future coverage uses the fitted gamma–Poisson distribution as a prior and conditions on the measured full-effort depth histogram. For a currently observed depth `d`, the accessible-site posterior rate is gamma with shape `α+d` and rate `α/μ+1`; additional depth at effort `t` is posterior-predictive negative binomial over exposure `t−1`. At zero observed depth, the prediction is also weighted by the posterior probability that the site belongs to the accessible component. Because no additional depth is possible at `t=1`, every forecast begins exactly at its measured breadth. One conditional distribution produces thresholds 1× through 100× jointly, guaranteeing that breadth cannot increase when the required depth increases and cannot decrease as sequencing effort increases.

Use the conditional model for exploratory scenarios through 100× effort. The grid is 0.01× to 2× and 0.25× thereafter. No coverage gives no model. The report displays breadth at the selected minimum depth without classifying references as complete or estimating required effort.

For every configured fraction `f` where `2f` is also present and at least four training fractions are available, a separate retrospective model treats the observed depth histogram at `f` as current and conditionally predicts `2f`. The report records maximum absolute breadth error across thresholds. This directly tests short-range extrapolation against reads withheld by deterministic hashing.

The model describes future reads from unchanged composition, library and bias. It cannot resolve true absence, divergence, mapping failure, primer dropout or extreme coverage heterogeneity from shallow sampling alone. Amplicon sequencing, major changes in library preparation and strongly structured coverage can violate the model. Extrapolations are not a substitute for empirical resequencing validation.

## References

* Chao et al. (2014), [Rarefaction and extrapolation with Hill numbers](https://doi.org/10.1890/13-0133.1): incidence-based richness estimation and extrapolation framework.
* Rodriguez-R et al. (2018), [Nonpareil 3](https://doi.org/10.1128/mSystems.00039-18): abundance-weighted metagenomic coverage, a different estimand from raw k-mer richness.
