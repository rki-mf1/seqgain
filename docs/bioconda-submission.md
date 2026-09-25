# Bioconda submission for SeqGain 0.3.0

The upstream release and the Bioconda recipe are separate: a Git tag supplies the source, while a pull request to [`bioconda-recipes`](https://github.com/bioconda/bioconda-recipes) publishes the package after review. The submission-ready recipe is at [`bioconda-recipe/meta.yaml`](../bioconda-recipe/meta.yaml).

1. Keep the public `v0.3.0` tag at its release commit. Do not move it after submission.
2. The recipe checksum was calculated from the public archive at `https://github.com/rki-mf1/seqgain/-/archive/v0.3.0/seqgain-v0.3.0.tar.gz`. To verify it, run:

   ```bash
   curl --fail --location --output seqgain-v0.3.0.tar.gz \
     https://github.com/rki-mf1/seqgain/-/archive/v0.3.0/seqgain-v0.3.0.tar.gz
   shasum -a 256 seqgain-v0.3.0.tar.gz
   ```

   Expected SHA-256: `fb8f3167f3cf33060ef878588a55839a308f81124510e71ffa5d03aa605f5c49`. Check that the download is the intended source archive, not an error page.
3. Sync your GitHub fork of `bioconda/bioconda-recipes`, create a branch from its current default branch, and copy `bioconda-recipe/meta.yaml` to `recipes/seqgain/meta.yaml` in that fork. The recipe maintainer is `simonhtausch`.
4. Run Bioconda lint/build tests if available, then open a pull request. Bioconda's CI will build and test the recipe; address any review or platform-specific failures there. The upstream GitLab CI checks the application but does not replace Bioconda's recipe checks.
5. After the PR is merged and the package is available, test a fresh installation on Bioconda, then update the README to remove its pre-publication warning.

Bioconda expects a stable source URL with a checksum, an appropriate license file, a build number of 0 for a new version, `noarch: python` for pure Python packages, command-line tests, and a `run_exports` pin. For a `0.x` version, Bioconda recommends a minor-version pin (`max_pin="x.x"`). The recipe follows those conventions.
