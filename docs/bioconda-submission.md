# Bioconda submission for SeqGain

The upstream release and the Bioconda recipe are separate: an immutable Git tag supplies the source, while a pull request to [`bioconda-recipes`](https://github.com/bioconda/bioconda-recipes) submits the package for review. Publication happens only after Bioconda accepts and builds the recipe.

1. Run the application tests and Conda package checks, then publish the release tag. SeqGain's first GitHub release is `v0.3.1`; the historical `v0.3.0` tag belongs to the former SeqScope project and must not be reused.
2. Download the tagged GitHub archive and calculate its checksum:

   ```bash
   curl --fail --location --output seqgain-v0.3.1.tar.gz \
     https://github.com/rki-mf1/seqgain/archive/refs/tags/v0.3.1.tar.gz
   shasum -a 256 seqgain-v0.3.1.tar.gz
   ```

3. Set that SHA-256 in [`bioconda-recipe/meta.yaml`](../bioconda-recipe/meta.yaml). This checksum update necessarily follows the source tag; do not move the tag to include its own archive checksum.
4. Create a branch from the current default branch of `bioconda/bioconda-recipes` and copy the completed recipe to `recipes/seqgain/meta.yaml`. The recipe maintainer is `simonhtausch`.
5. Run the recipe build and Bioconda lint checks, then open a pull request. GitHub Actions checks the upstream application; Bioconda's CI additionally checks its packaging requirements and supported build platforms.
6. After acceptance and package publication, test a fresh Bioconda installation before removing the README's pre-publication warning.

The recipe uses a tagged source archive, SHA-256 verification, build number 0, `noarch: python`, the GPL-3.0-only license, import and CLI checks, and a minor-version `run_exports` pin for this `0.x` package. See the [Bioconda recipe guidelines](https://bioconda.github.io/contributor/guidelines.html).
