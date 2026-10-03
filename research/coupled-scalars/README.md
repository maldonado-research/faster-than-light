# Coupled preferred-time scalar verification

2 October 2026 Pacific (3 October UTC). These are registered internal synthetic mathematical checks of an assumed scalar model. No physical channel, detector, observation, novel theory or external peer review is claimed. See the [derivation and scope](../../docs/COUPLED_SCALAR_RESPONSE.md).

The stable mass region is `g^2<=m^2*M^2` on R^d or periodic domains admitting k=0. Stable massive examples have a rigorously positive **1+1-dimensional** chi-source/chi-response in the intercone at sufficiently small preferred times. The same-field massless convolution coefficient starts at order g^2 and vanishes quartically near the faster front; its 3D radial coefficient has cubic onset. A massless nonzero-g full system is unstable. No 3D positive-response inequality is obtained by differentiating the 1D bound.

## Registered results

| Route | Evidence and final checks |
| --- | --- |
| Primary | [Registration](REGISTRATION.md), preserved [original](REGISTRATION_ORIGINAL.md), [implementation protocol](root/IMPLEMENTATION_REGISTRATION.md), [results](root/RESULTS.json) and [source binding](root/RESULTS_BINDING.json). 2,450 case evaluations: 2,000 mode/energy, 200 exact transform, 120 piecewise support/front, 60 exact spatial-integral, 20 high-k residue, 16 mass boundaries, 24 coincident cones, six k=0 transforms and four parameterized pathological cases; plus 960 interior grid and 16 scaling-ratio diagnostics. |
| Separate implementation | [Registration and amendment](independent/REGISTRATION.md), [final results](independent/RESULTS.json), [audit](independent/AUDIT.md) and [manifest](independent/MANIFEST.json). 10,415 constraints, including 54 stable full-response certificates, spectral versus matrix-series retarded checks, exact normalization/3D radial checks and low-energy elimination remainders. The healthy 63-constraint gate passes; four actual injected defects fail it. |
| Massive kernel | [Registration and amendments](kernel/REGISTRATION.md), [results](kernel/RESULTS.json), [derivation](kernel/DERIVATION.md) and [audit](kernel/AUDIT.md). 36 exact triangle moments, 48 exact transform identities, 14 configured massive-cell evaluations covering 12 distinct cells, 28 grids/56 scalar integrals, four massless normalizations, 14 exact full-response bound cells, eight monotonicity/eight front-power comparisons, two slope checks and three stability controls. |

All final routes passed their registered criteria. Counts refer to different kinds of synthetic constraints and are not numbers of experiments. The primary maximum normalized eigen residual is 3.98e-16; its raw characteristic residual scales with the squared matrix norm and is reported separately. Massive quadrature at 64/128 Simpson subdivisions changes by at most 3.06e-11 relatively, within the declared 1e-9 tolerance. Numerical quadrature is a check of the order-g^2 coefficient; the full-response bounds come from the convergent series proof.

In the normalization of the 1D impulse kernel, the stable example `a=2,m=M=1,g=1/2,t=1/10,r=3/20` has order-g^2 response approximately `1.808015842e-9`. The exact full response is bounded approximately between `1.806000133e-9` and `1.808449077e-9` by the documented rational inequalities. Those small numbers have no assigned laboratory units, signal-to-noise or apparatus-energy interpretation.

## Reproduce

Python>=3.10 and the standard library only. From the repository root, choose fresh output paths; existing files are refused:

```bash
python3 research/coupled-scalars/root/reproduce.py --out /tmp/ftl-coupled-root.json
python3 research/coupled-scalars/independent/check.py --out /tmp/ftl-coupled-independent.json
python3 research/coupled-scalars/kernel/reproduce.py --out /tmp/ftl-coupled-kernel.json
```

The preserved source-bound result files are not overwritten by these commands. UTC run times differ on reproduction. Input parameters, seeds, tolerances, source/registration hashes, counts, failures and numerical extrema are stored in the JSON files. The original central-registration hash binds the separately preserved original text; the amended registration is also hashed in the final runs.

## Disclosed corrections

The initial separate implementation passed 10,414 checks but its four claimed mutation rejections were only comparisons with incorrect values, and two pathology checks were literal sign calculations. That certification was insufficient and is withdrawn. The amended implementation passes the same parameterized checker over the healthy candidate and four actual injected candidates, uses shared pathology formulas, and adds one healthy gate for final 10,415. [Initial source](independent/INITIAL_CHECK.py), [results](independent/INITIAL_RESULTS.json), registration/audit snapshots and final amendment remain available. The initial result is historical evidence, not the strengthened certification.

The kernel's first source restricted output to its own work directory. A source-only portability amendment now permits an explicit fresh external output and refuses overwriting source/canonical inputs. Initial source/results are preserved; the complete mathematical payload agrees after the source-only rerun. The pre-run description incorrectly called 14 evaluations 13 distinct cells; the disclosed correction is12 distinct, with no numerical count or tolerance change.

Analytical review clarified internal energy versus the forced canonical Hamiltonian, the domain assumed by the k=0 stability criterion, and the canonically normalized derivative-truncated mass. The [skeptical scope audit](review/SKEPTICAL_SCOPE_AUDIT.md) and final source binding document these boundaries. These are internal reviews and disclosed corrections, not external scientific replication.

The proposed Zenodo v1.2.2 ZIP remains frozen to preceding methods and contains none of this round. Missing private v6.3 code/proofs and strict observational inputs remain unresolved. A finite-energy smooth source/response calculation, real detector identification, noise/energy budgets, transferred Lorentz-violation constraints and quantum/UV/gravity completion remain future work.
