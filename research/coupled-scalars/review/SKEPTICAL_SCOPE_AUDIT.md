# Coupled-scalar analytical and claim-boundary audit

Reviewer: workflow_review. Initial review: 2026-10-03 UTC / 2026-10-02 Pacific. Read-only inspection of `COUPLED_SCALAR_RESPONSE_DRAFT.md`, `REGISTRATION.md` and repository guidance; no root implementation read or executed. These are independent analytical checks, not observations, external peer review or a literature novelty assessment. This review does not alter the frozen proposed v1.2.2 package or authorize a publication action.

## Findings requiring precise wording

1. **Mass coefficient versus canonical mass.** Eliminating phi gives `D_s-g²(D_f^R)^-1`. The displayed heavy-mass expansion is correct: `Z_t=1+g²/m⁴`, `Z_x=1+a²g²/m⁴`, and potential coefficient `C=M²-g²/m²`. Calling C the effective mass squared without qualification is imprecise. After normalizing chi's kinetic term, the truncated low-energy mass squared is `C/Z_t`, and the effective squared propagation coefficient is `Z_x/Z_t`. Higher-derivative corrections remain. Prefer “unnormalized mass coefficient” for C or show the division by Z_t explicitly.
2. **The Hamiltonian iff statement needs a spatial-domain assumption.** Mass-matrix positive semidefiniteness is sufficient on every usual domain. It is necessary for the field Hamiltonian on R^d when arbitrarily broad compactly supported configurations are allowed, or on periodic boxes admitting constant k=0 configurations. A fixed bounded Dirichlet domain can have a stabilizing gradient/Poincare gap, so a negative mass direction need not make its restricted energy unbounded below. State the unrestricted/low-k setting before claiming the precise iff criterion and low-k tachyon control.
3. **Energy of an externally forced model.** The balance `partial_t e+div S=J chi_t` is correct for internal field energy e, whose source-free integral is the preferred Hamiltonian. The canonical Hamiltonian of the written L including prescribed `+J chi` is `e-J chi`; its balance instead involves `-J_t chi`. Naming e “internal field energy density / source-free Hamiltonian density” avoids suggesting the externally driven system has conserved or bounded total source-plus-apparatus energy.

## Independently checked derivations

The equations of motion, sign of the mixing potential, preferred flux, mass determinant criterion, mode eigenvalues and trace/product are correct. Subtracting k² times identity from the dispersion matrix leaves the stable mass matrix plus `diag((a²-1)k²,0)`, so `omega_±²>=k²` follows under the assumed unrestricted stability criterion. Zero-frequency saturation is a neutral boundary; it should not be advertised as strictly positive frequency or uniformly bounded homogeneous amplitude. A zero mode can have a linear-in-time homogeneous solution.

The exact retarded elimination has the sign

`phi=-g G_f*chi`,
`chi=G_s*J+g² G_s*G_f*chi`.

Thus the cross response is odd in g and starts at `-g G_f*G_s`, whereas the same-field response is even and has the stated Volterra series. Its intercone g² coefficient is a Born/Duhamel coefficient, not by itself the full stable response.

An independent Laplace/Fourier derivation confirms the H and K formulas. Write `S=s²`, `d=a²-1`, `Q_f=S+a²k²`, `Q_s=S+k²`. Then

`1/(Q_f Q_s) = a²/(d S Q_f)-1/(d S Q_s)`,

`1/(Q_f Q_s²) = a⁴/(d² S² Q_f)-a²/(d² S² Q_s)-1/(d S Q_s²)`.

Inverse time integration of the free kernels gives precisely the draft's q_a/q_1 expressions. Integrating a massless free kernel over x gives t independently of c, so spatially integrated H and K are respectively `t³/3!` and `t⁵/5!`. Coincident-speed kernels `theta(t-r)(t²-r²)/8` and `theta(t-r)(t²-r²)²/128` have the correct normalization and must be handled separately from singular distinct-speed denominators.

## A complete one-dimensional convergence and positivity argument

For every free 1+1 massive kernel, `|J0(z)|<=1` gives the spatial bounds `||G_c(t,.)||_infinity<=1/(2c)<=1/2` and `||G_c(t,.)||_1<=t`, for c>=1. In a causal convolution of 2n+1 kernels, retain one factor's infinity bound and use the other 2n spatial L1 bounds. The remaining time-simplex integral is

`(1/2) integral_{tau_i>=0,sum(tau_i)<=t} product(tau_i) d^{2n}tau = t^(4n)/(2(4n)!)`.

After multiplying by `|g|^(2n)`, the bound is summable and uniform on every bounded time interval. The resulting causal series solves the Volterra integral equation; uniqueness follows from the same iteration bound applied to a homogeneous difference. Therefore the series is the full retarded response, not only a formal perturbation expansion. Its absolutely convergent kernel representation can be used for the subsequent lower bound. Global absolute magnitude bounds do not imply stable dynamics for parameter values outside the mass positivity region.

For `m t<=1` and `M t<=1`, all segment durations are <=t and the respective J0 arguments are <=1. Each free kernel is nonnegative, with lower factors `1-m²t²/4` or `1-M²t²/4`. Every even-mixing convolution is therefore nonnegative. In `t<r<a t`, the n=0 ordinary kernel vanishes, but the n=1 term is bounded below by

`g² (1-m²t²/4)(1-M²t²/4)² (a t-r)^4/[48 a(a²-1)²] >0`

when g is nonzero. This proves a full stable same-field mathematical response in the small-time intercone. The example `a=2,m=M=1,g=1/2` satisfies strict mass positivity. Quadrature of the massive n=1 term can verify that coefficient, but should not be labelled numerical summation of the full response; the full-response inference is analytical.

## Dimension, front and apparatus boundaries

The dimensional derivative relation gives the correct massless radial intercone coefficients H3 and K3. The same-field radial coefficient has cubic onset, while K1 has quartic onset. Derivatives must be distributional where fronts/contact terms occur.

**Keep the above positive-kernel proof and pointwise lower bound explicitly 1+1 dimensional.** In 3+1 dimensions, differentiating the free massive 1D kernel produces a delta front and a negative interior J1 tail:

`G_c,mu^(3)=delta(t-r/c)/(4 pi c² r)-theta(t-r/c) mu J1(mu tau)/(4 pi c³ tau)`,

where `tau=sqrt(t²-r²/c²)` and r>0. Consequently termwise positive Volterra kernels do not transfer automatically to 3D. Nor may one differentiate a scalar inequality to obtain a derivative inequality. A complete 3D front argument may instead use the exact dimensional relation, continuity/support at the outer cone, or a separately controlled near-front expansion. The draft currently calls H3/K3 perturbative radial coefficients, which is the correct limitation.

The outer causal-support bound follows from the local hyperbolic system and causal convolution supports. For 1D, the positive response arbitrarily near the fastest cone establishes a nonzero mathematical front at a for the complete continuum operator. Finite bandwidth and a low-energy derivative expansion do not fix a physical ultraviolet front. Nonlocal spectral projections can have separate spatial tails; only the complete local matrix response has the stated compact causal support.

Common future-directed preferred time excludes a finite closed signaling/relay loop only when every channel and relay respects the same global ordering. The scalar coupling does not restore the earlier reciprocal sender-frame rule or permit reuse of the original J/occupancy-to-speed diagnostic without a new population/detector model. No PTA/ringdown calibration follows.

The slower-field slope and induced Lorentz-violating derivative coefficients are correct under the stated assumptions. They also mean ordinary matter cannot simply be identified with chi while leaving rods, clocks, source control and detector timing unchanged. The small-k slope is not the k=0 group velocity of a gapped massive mode. Positive source-free energy is separate from physical apparatus work, pulse cost, backreaction, noise and detection thresholds. The draft's prescribed scalar forcing has no derived photon/matter identification or measurement observable.

## Claim/reporting status

The registration openly discloses prior analytical knowledge, internal reviewers and its pre-run extension. The draft correctly excludes physical FTL detection, external peer review, blind discovery preregistration, novelty, strict empirical v1.3.0 and recovered v6.3/observational inputs. No numerical pass counts appear prematurely in the examined draft. Final source-bound results, amendments and failed/skipped outcomes remain to be checked when supplied.

The owner-facing Zenodo guide is reviewed separately; new coupled-scalar work must remain outside the frozen proposed v1.2.2 ZIP/JSON. Valid complete existing-draft owner UI publication is the owner's decision and is not implied by this analytical audit.

## Final local-results follow-up — 3 October UTC

The current draft resolves all three initial wording findings: it specifies unbounded/periodic low-k domains and the Dirichlet exception, distinguishes internal field energy from `e-J chi`, and displays the canonically normalized heavy-field mass. It also explicitly limits the massive positive-response inequality to 1+1 dimensions and warns that a 3+1 negative Bessel tail prevents differentiation of that inequality into a positive-response bound.

The final saved reporting/binding review passes 24/24 artifact and claim checks, recorded separately in `FINAL_REPORTING_AUDIT.json`. No scientific suite was reexecuted. The root's 2,450 cases have empty failures, with 960 internal grid comparisons and 16 residue-ratio comparisons reported separately. Its raw characteristic residual exceeds an absolute 1e-10; the disclosed pre-run normalized residual is below the declared threshold, and the run report does not represent the raw residual as an absolute pass.

Kernel verification reports 14 configured massive/front evaluations, **12 distinct parameter cells**, 28 grids and 56 scalar integrals. Two front points repeat stable-cell points. The original registration's descriptive 13-distinct count is preserved and corrected by a disclosed amendment; no scientific case or threshold was removed. The output-portability source change and final rerun are also disclosed. Comparing initial and final JSON after excluding run/source/registration hashes shows identical mathematical payloads. Two-resolution agreement is numerical convergence evidence, not a rigorous Simpson-error certificate; the stable full-response inference rests on the analytical bounds.

The independent final result records **10,415 checks**, an empty failure list and 54 strictly positive exact full-response lower-bound certificates. It preserves the initial 10,414-check artifacts and openly withdraws their weak purported mutation-rejection/control certification. A disclosed post-run amendment adds genuine candidate substitutions and one healthy-candidate gate, while strengthening the shared negative controls. The final result records all 63 healthy inner constraints passing and four injected defects rejected; these inner gate evaluations are not inflated into the headline total. This reviewer checked result/registration/manifest consistency and file hashes, not executable gate semantics; the assigned implementation reviewers own that separate check.

Root, kernel and independent final script hashes match their saved results, and the independent current/initial manifest matches the actual artifact hashes. Only script bytes were hashed; the root implementation contents were not inspected. Final public-source copying, relative links and document pass-count wording were still pending when this follow-up was recorded and must be bound after staging. No physical FTL observation, apparatus feasibility, novel theory or external peer review follows from these internal conditional checks.
