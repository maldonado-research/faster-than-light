# Reproduce the recovered public synthetic smoke path

This is a current software check of unchanged public v6.0 code, using explicitly synthetic fixtures. It does not reproduce missing v6.3 guarantees or any physical observations. [Live archive verification](../../docs/LIVE_ARCHIVE_VERIFICATION_20261001.md) records official sources and checksums; [VALIDATION.json](VALIDATION.json) records current output validation and actual library versions.

Obtain FTL_A2_STRICT_PROTOCOL_FREEZE_KIT_v6_0.zip from [record 20218030](https://zenodo.org/records/20218030). Verify its 1,502,150-byte length and declared MD5 d9fb0abc2ae08c06bf2091171273191f before use, and verify safe ZIP paths before extracting outside your Git checkout. Change into the extracted kit directory containing SCRIPTS/, VENDOR/ and EXAMPLES/.

The current run used Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, Matplotlib 3.10.8, SciPy 1.17.0 and scikit-learn 1.8.0. The published archive has no locked dependency specification. This records the tested environment rather than promising arbitrary versions work.

From that kit directory, select a new output directory and run:

```bash
set -e
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/ftl-a2-v6-cache/matplotlib
export XDG_CACHE_HOME=/tmp/ftl-a2-v6-cache
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME/fontconfig"
timeout 180s python3 SCRIPTS/ftl_a2_v6_strict_release_runner.py \
  --provider-bundle EXAMPLES/FTL_A2_CROSSFIT_SYNTHETIC_PROVIDER_BUNDLE_v5_9.zip \
  --metadata EXAMPLES/provider_metadata_synthetic_v5_9.json \
  --outdir /tmp/ftl-a2-v6-synthetic-smoke \
  --label ftl_resumption_synthetic_v6 --allow-synthetic \
  --bootstrap 150 --pci-bootstrap 12 --n-per-class 800 \
  --folds 3 --seed 60 --wcap 10
```

The wrapper exit status alone is insufficient: it does not propagate failures of its children. Inspect the new FTL_A2_V6_PROTOCOL_RUN_INDEX.json and each run's logs/reports. Require all three child exit codes to be zero, provider pass, COMPLETE reports, 150 finite radicand draws and 12 finite PCI bootstrap draws per configuration. Recompute phase masses and thresholds from current draws and report invalid radicand mass separately. Do not use archived output timestamps as current-run evidence.

The current audited run passed these checks, with zero invalid radicand draws. Its three classifier stability classifications were STABLE_5PCT. Stability is not an accuracy certificate; later source-reported classifier and provenance defects remain relevant. A failure or changed status in a new run must be reported directly.
