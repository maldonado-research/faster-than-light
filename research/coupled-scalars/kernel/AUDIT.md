# Independent bounded verification of the 1+1 dimensional kernel

This is an internal implementation and mathematical audit of the assumed coupled-scalar model, not an external replication, physical experiment, novelty assessment, or detector calculation. The algebra in DERIVATION.md was known before registration. REGISTRATION.md discloses that order and defines the cases, counts and acceptance criteria. The root registration is in ../REGISTRATION.md. This implementation and its results are separate from the root and independent-review implementations.

## Outcome and reproduction

All configured mathematical checks passed with no recorded failures. Python 3.12.14 and the standard library were used. Run time was approximately one second. The final inspected output was first written to `/tmp/ftl-kernel-final.json`, then explicitly copied to RESULTS.json. Its recorded UTC time is `2026-10-03T00:57:37.306010+00:00`.

Keep the containing registration beside the root, kernel and independent directories. From any working directory, run the kernel script with an explicit fresh writable output path, for example:

```sh
python research/coupled-scalars/kernel/reproduce.py --out /tmp/ftl-kernel-check.json
```

The path shown assumes the public repository layout. An existing output, symlink, canonical result, source or registration path is refused; output creation uses exclusive mode. A new file is required for each reproduction. The command reads both its own registration and its parent's registration to bind their hashes. A successful run exits zero; an exception or failed assertion is recorded and exits nonzero.

## Configured checks

| Check group | Completed count |
| --- | ---: |
| Exact triangle moments, p,q=0,...,5 | 36 |
| Exact transform cases | 24 |
| Exact H/K partial-fraction identities | 48 |
| Massive/front evaluations | 14 |
| Massive/front Simpson grids | 28 |
| Massive/front scalar integrals | 56 |
| Massless normalization grids | 2 |
| Massless normalization checks | 4 |
| Exact-rational full-response-bound cells | 14 |
| Adjacent normalized-ratio monotonicity checks | 8 |
| Near-front power checks | 8 |
| Final front-slope checks | 2 |
| Stability controls | 3 |

Counts in this table include cases, grids and values and must not be summed as independent assertions. The 14 massive/front evaluations cover 12 distinct cells. The stable grid has `a=2,m=M=1,g=1/2`, `t={.01,.05,.1}` and `r/t={1.25,1.5,1.75}`. The five additional front evaluations have `t=.1` and `Delta/t={.5,.25,.125,.0625,.03125}`. Their first two points repeat two stable-grid cells.

The massive triangle calculation used composite Simpson resolutions 64 and 128. The maximum relative change was `3.056596934885074e-11`, below the registered `1e-9` criterion. Each normalized coefficient also met its analytic mass bound with the registered `2e-12` numerical slack. Agreement of two quadrature resolutions is a numerical consistency check, not a rigorous quadrature-error enclosure.

The Bessel factors were evaluated by twelve-term entire alternating series with squared argument at most one. The analytic first-omitted-term bounds are `1/[4^12*(12!)^2]` for J0 and `1/[4^12*12!*13!]` for `2J1/z`. Those bounds apply to series truncation, separately from floating roundoff and quadrature error. The exact massless normalization values are the triangle moments 1/2 and 1/24; the two configured resolutions met the `1e-13` absolute criterion.

For the last front point, the measured normalized H slope was `0.008300700631167501`, versus `m^2 t/(6a)=0.008333333333333333`, with relative difference `0.003915924259899825`. The K slope was `0.004986946783489543`, versus `m^2 t/(10a)=0.005`, with relative difference `0.0026106433020913963`. Both meet the registered one-percent criterion. These finite-gap checks support the stated asymptotic coefficients; the coefficients themselves follow analytically from triangle moments.

At `a=2,m=M=1,g=.5,t=.1,r=.15`, the computed leading coefficient `g^2 K` is approximately `1.808015842024381e-9`. The rigorous analytic full-response interval, independently of numerical quadrature, is

\[
\frac{9587}{5308416000000}\le G_{\chi\chi}\le
\frac{38707200041}{21403533308129280000},
\]

approximately `[1.8060001326195988e-9, 1.8084490763166945e-9]`. The lower bound uses the massive leading-kernel inequality and nonnegativity of every even Volterra term at these short times. The upper bound uses the massless majorant and convergent higher-mixing bound. Thus the interval concerns the full assumed kernel; the reported leading numerical coefficient is not represented as a numerical solution of that full kernel.

The stable control has potential determinant 3/4 and minimum squared frequency 1/2 at k=0. Both nonzero-mixing massless controls, g=+/-1/2, have minimum squared frequency -1/2 and are classified as unstable. No unstable full model is used as the stable-response example.

## Source and result binding

| Record | SHA-256 |
| --- | --- |
| Final reproduce.py | `4bdd99418bdd372703f03ea275b7f2d840ebeeb61469ae6e1403c323f1d3d5f4` |
| Amended kernel REGISTRATION.md | `9cb442a3ba3a6972d10eda6bf08a38a1fe13f6b400999f14f5b3327507a1c58a` |
| Root ../REGISTRATION.md at final run | `23c5d75de7f38c9f97388fcc4fb941e964449766923164a0a10fc7c55c2fc841` |
| Final RESULTS.json | `42db6eede01ec058dd2ba97c443fe014799483943b9bc090139a44ed515d5eb1` |
| INITIAL_REPRODUCE.py | `1830892f8b9e3c0c6b1a449ba5224cfe34459e4e60e049931bfa148dbe194388` |
| INITIAL_RESULTS.json | `d4dd2420ef4b47260281858f158b8e4aea4f09a1e6585439b89fbc151fecfb3c` |
| Kernel registration before amendment | `59d2de44809f56a0e8530cc7cacfe30877b75a867fc7d2a3dc752f5cff0fbee5` |

The initial script and results were preserved byte-for-byte before the source-only amendment. The first successful run was at `2026-10-03T00:50:18.368314+00:00`. Its directory-only output restriction was inherited from the development task and prevented a portable `/tmp` reproduction. After inspection, root authorized a guard correction. The correction permits an explicit fresh writable output while refusing replacement of existing files and protected records; it changes no formula, mathematical count, parameter, tolerance or resolution.

The original registration's claim of 13 distinct cells was a descriptive count error: the configured evaluations actually cover 12. The original text is preserved, and the amendment identifies both repeated cells. This correction was made before the final run and changes no configured evaluation. The final registration hash binds both the portability amendment and count correction.

The final rerun's complete mathematical payload was compared with the initial output. It is exactly unchanged. Only `run_utc`, `script_sha256` and `registration_sha256` differ. The root-registration hash remained unchanged. Two separate post-run output-guard checks confirmed refusal to replace the final temporary output or the script itself, with each target byte-identical afterwards; they did not execute new mathematical checks and are not included in the registered mathematical counts.

## Scope and limitations

The exact support and positivity claims come from continuum retarded convolutions and a convergent Volterra series. The numerical calculation uses a bounded triangle, not a cutoff Fourier reconstruction, and cannot create an instantaneous-support claim. The table verifies 1+1 dimensional kernels only; it does not establish a 3+1 dimensional numerical result. Any dimensional extension requires its own derivation and verification.

The fields, preferred frame, fast characteristic speed and local mixing are assumptions. Chi is a field label, not an identified real material or calibrated detector observable. The impulse forcing supplies no finite source construction or energy budget. Stability of this classical quadratic model is distinct from a microscopic or ultraviolet completion. There is no experimental evidence, usable communication device, measured signal size, novelty claim, or external peer-review endorsement in this audit. No Git, website or Zenodo publication was performed as part of these runs.
