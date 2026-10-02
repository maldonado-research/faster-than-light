# Finite proper reply delay: replication and sensitivity protocol

Recorded 2026-10-01 Pacific (2026-10-02 UTC). This is a subsequent mathematical methods exercise, separate from the September 8 research snapshot and September 30 documentation edition.

Prior knowledge: the immediate-reply formula was previously reproduced. The delay formulas were derived algebraically and checked by an internal agent before this protocol was written. This is a disclosed replication/sensitivity registration, not a blind discovery preregistration. No novelty or observational claim is proposed.

Model: Minkowski spacetime; collinear inertial Alice and Bob; coincident origins; Alice sends at T>0; equal hypothetical controllable speed Wc in each sender's frame; Bob waits proper time tau>=0 before replying. Define beta=v/c and d=tau/T. The severity population uses uniform beta on [0,1) and a fixed d, independent of W. This is an explicit model extension.

Required checks, fixed before executing the accompanying implementation:

1. Compare 5,000 deterministic direct Lorentz-event calculations (seed 20261001) with the delay formula; relative/absolute tolerance 1e-10. Add exact Fraction comparisons for rational beta and complementary square root pairs.
2. Evaluate 30 thresholds using 70-digit Decimal arithmetic (W=1.001,1.1,2,10,100; d=0,0.1,0.5,1,10,100). Require residual below 1e-45 and correct strict signs immediately on both sides.
3. Verify the exact example W=2, beta=12/13, d=1/2 gives t2/T=85/98 and maximum d=24/35.
4. Calculate uniform-beta infinite-W ceilings for d=0,0.1,0.5,1,2,10,100,1000. Compare Simpson resolutions 4096 and 8192 at relative tolerance 1e-9. Require positivity, strict decrease, d=0 recovering 1/5 within 1e-12, and d=1000 giving d^2 J within 1e-6 of 1/12.
5. Compare finite-W severity at W=1.1,2,10,100,1e6 for d=0,0.1,0.5,1,2,10. Require monotonic increase with W, values below the corresponding ceiling, and the last value within relative 1e-4 of that ceiling.
6. A separate internal audit checks derivation, integration measure, parameter assumptions, thresholds and claim boundaries. A zero test count, skipped run or failed check does not establish success. Failed checks preserve a nonzero exit status.

All inputs are deterministic mathematical parameters. There are no posterior observations, private raw data or synthetic observational fixtures in this calculation. Reproduce from this directory with `python3 reproduce.py --out /tmp/ftl-response-delay-results.json`, then inspect that current result. Do not overwrite archived source equations or treat changes to this conditional severity as physical evidence.
