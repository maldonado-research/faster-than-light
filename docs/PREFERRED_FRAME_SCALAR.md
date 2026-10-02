# A fixed-preferred-time scalar propagation exercise

Ricardo Maldonado · subsequent conditional-theory exercise · 2026-10-01 Pacific (2026-10-02 UTC)

The earlier antitelephone rounds assumed a controllable speed in each sender's frame. This note instead specifies an illustrative free action with one preferred time. It passes the stated classical free-field energy and stability checks and supplies conditional timing predictions. It is a familiar class of toy model, not a new discovery, Ricardo's established hidden-sector theory, a physically realized FTL channel, an observational result or external peer review. The historical snapshot remains unchanged.

## Explicit model change

Use c=1 and Minkowski metric eta=(-,+,+,+). Choose a fixed, globally constant unit timelike background u, and h^{ab}=eta^{ab}+u^a*u^b. For a real scalar phi, assume

$$
\mathcal L=\frac12(u^a\partial_a\phi)^2
-\frac{a^2}{2}h^{ab}\partial_a\phi\partial_b\phi
-\frac{m^2}{2}\phi^2+J\phi,\qquad a>1,\quad m^2\ge0.
$$

The symbol a is the dimensionless preferred propagation speed. J is prescribed mathematical forcing, with no emitter/detector or interaction action derived. Ordinary rods, clocks and relays obey the Minkowski metric by assumption. u is not dynamical, and no gravitational or quantum completion is specified.

Let P^{ab}=u^a*u^b-a^2*h^{ab}. Variation gives partial_a(P^{ab}partial_b phi)+m^2*phi=J. In preferred coordinates u=(1,0,0,0),

$$
\boxed{\partial_t^2\phi-a^2\nabla^2\phi+m^2\phi=J.}
$$

One can transform coordinates and u together. A chosen fixed-u state nevertheless breaks boost symmetry; coordinate covariance does not make the physical propagation rule isotropic in every observer's rest frame.

## Preferred energy and stability

With J=0, the canonical momentum is pi=partial_t phi, and

$$
e=\frac12[\pi^2+a^2|\nabla\phi|^2+m^2\phi^2]\ge0,
\qquad \mathbf S=-a^2\pi\nabla\phi,
\qquad \partial_t e+\nabla\cdot\mathbf S=0.
$$

For sufficient decay or appropriate boundary conditions, H=integral(e d^3x) is conserved and bounded below. Nonnegative includes the massless constant zero mode where allowed by the boundary conditions. With forcing, partial_t e+div S=J*partial_t phi; source power and backreaction are not supplied. Positivity of the free Hamiltonian is not positivity of an arbitrary prescribed-source or interacting system.

Plane waves obey

$$
\omega^2=a^2|\mathbf k|^2+m^2,\qquad
v_{\rm phase}=\frac\omega{|\mathbf k|},\qquad
\mathbf v_{\rm group}=\frac{a^2\mathbf k}{\omega}.
$$

All preferred free frequencies are real. Phase speed is at least a, group speed is at most a, and the massless speeds both equal a. For m>0 the group speed exceeds ordinary light speed precisely when |k|^2>m^2/[a^2*(a^2-1)]. A negative m^2 produces a low-k instability; a negative gradient coefficient produces a high-k instability; reversing the overall action sign leaves the source-free equation but makes the Hamiltonian unbounded below. Those controls are detected separately in the reproduction.

For a nonzero preferred positive-frequency mode,

$$
\eta^{ab}k_a k_b=-\omega^2+|\mathbf k|^2
=-(a^2-1)|\mathbf k|^2-m^2<0.
$$

Its wave covector is Minkowski timelike despite potentially superluminal group speed. Treating every faster-cone excitation as a Lorentz-invariant spacelike-momentum tachyon would misidentify this changed model. The JSON's positive `boost_norm` quantity is minus the displayed eta contraction.

## Retarded front and its ultraviolet limit

The principal cone is |dx|<=a*dt, oriented by dt>0. The energy bound |S·n|<=a*e implies zero energy remains outside an initially empty expanding region with boundary speed a. Thus compact initial data and retarded forcing have no response outside their future preferred cones.

For the complete 3+1-dimensional massless continuum PDE, its retarded fundamental solution at r>0 is

$$
G_{\rm ret}(t,r)=\frac{\delta(t-r/a)}{4\pi a^2r}.
$$

The front is therefore a in this explicit complete free model. A mass term changes the interior response, not the characteristic cone. A compact 1D massless wave phi=F(x-a*t), F(z)=(1-z^2)^3 for |z|<=1 and zero outside, provides a finite-energy check: its energy is a^2*(1024/385) in the stated unit-amplitude, unit-width normalization, and S/e=a where nonzero. Two quadrature resolutions verify that energy; the domain-of-dependence conclusion also has the analytic proof above.

**A low-energy effective action alone does not determine the physical front.** The stated front conclusion assumes the operator remains valid at arbitrarily high frequency. An unknown cutoff/completion, a finite background region or a different operational coupling can change that conclusion. The delta-source solution is a distribution, not a finite-energy emitter or practical signal design. Scaling a classical field amplitude scales its energy quadratically; that supplies no minimum apparatus cost or production mechanism.

## Boosted initial surfaces and apparent negative coefficients

For t'=gamma(t-v*x), x'=gamma(x-v*t), |v|<1, the principal coefficients are

$$
A=P^{t't'}=\gamma^2(1-a^2v^2),\quad
B=P^{t'x'}=\gamma^2v(a^2-1),\quad
C=P^{x'x'}=\gamma^2(v^2-a^2),\quad AC-B^2=-a^2.
$$

The t'=constant hyperplanes are Cauchy surfaces for this preferred cone when |v|*a<1. Equality is characteristic; above it they are timelike relative to the propagation cone. The negative A on such a surface does not, by itself, identify a ghost in the preferred free theory.

Arbitrary t'-slice Fourier data would satisfy

$$
A\omega'^2-2Bk'_x\omega'+Ck_x'^2-a^2|k'_\perp|^2-m^2=0,
$$

with discriminant divided by four D=a^2*kx'^2+A*(a^2*|kperp'|^2+m^2). For a=2,v=3/4,kx'=0,|kperp'|=1,m=0, D=-80/7. Unbounded transverse wave numbers then produce unbounded imaginary-frequency growth in an arbitrary-data t' evolution problem. This is a 3+1-dimensional wrong-slicing obstruction. Such data are not physical finite-energy states specified on preferred Cauchy surfaces; the transverse argument does not assert universal ill-posedness for every massless 1+1-dimensional timelike-slice construction.

Physical preferred on-shell modes transform to omega'=gamma(omega-v*kx)>0 for every ordinary |v|<1, since omega>=a*|k|. Define spatial momentum P_x=-integral(pi*partial_x phi) on a preferred Cauchy slice. The observer translation charge is Q=gamma(H-v*P_x)>=0, with square completion

$$
e+v\pi\phi_x=
\frac12[(\pi+v\phi_x)^2+(a^2-v^2)\phi_x^2
+a^2|\nabla_\perp\phi|^2+m^2\phi^2]\ge0.
$$

This conserved charge differs from a formal canonical surface integral on a non-Cauchy t' hyperplane. Even a physical mode can have signed flux there. The exact example a=2,v=4/5 transforms (omega,kx)=(2,1) to (omega',kx')=(2,-1), with positive omega' but A*omega'-B*kx'=-2. The independent script verifies this distinction.

## Common time and the existing loop diagnostic

All operational preferred-retarded signals increase the same t. Ordinary future-directed matter/light relays also increase it. Nonnegative waiting cannot reverse that ordering, so a finite nontrivial chain cannot close or reach an earlier point on the sender's worldline. This requires every available channel to respect that common global ordering; it is not a theorem about unspecified interactions, curved spacetime, variable u or arbitrary shortcuts.

For Alice at rest in the preferred frame, Bob receding at v, an initial send at T>0 and Bob's proper wait tau=d*T, worldline intersections give

$$
\boxed{R_{\rm preferred}=\frac{t_2}{T}
=\frac{a+v}{a-v}+\gamma d\left(1+\frac va\right)\ge1.}
$$

At a=2,v=9/10,d=0, the [reciprocal sender-frame rule](FINITE_RESPONSE_DELAY.md) gives 76/121<1, whereas this preferred-frame rule gives 29/11>1. The two answers arise from different propagation assumptions.

If the positive loop strength is computed from this actual preferred prescription, P_+=max(1-R_preferred,0)=0. Its occupancy, first positive moment and squared severity J are all zero for any population supported on these assumptions. Consequently the proposed x=20*DeltaNeff*J channel cannot constrain a here. This is a diagnostic-scope limitation, not physical evidence for FTL or a result about the original sender-frame model. Other constraints and the missing physical/observational bridge remain relevant.

## Conditional laboratory signature

Assume a co-moving emitter/receiver with an ordinary Lorentz rest baseline L, preferred velocity v along x, massless controllable phi pulses and calibrated negligible reply latency. These apparatus and coupling assumptions have not been derived. The sending clock's proper roundtrip times for parallel and transverse baselines are

$$
\boxed{\tau_\parallel=\frac{2aL(1-v^2)}{a^2-v^2},\qquad
\tau_\perp=2L\sqrt{\frac{1-v^2}{a^2-v^2}}.}
$$

They recover 2L/a at rest and 2L for light-speed a=1. At small preferred velocity,

$$
\tau_\parallel-\tau_\perp
=-\frac La(1-a^{-2})v^2+O(v^4).
$$

This orientation-dependent roundtrip difference avoids relying solely on one-way clock synchronization. Restoring units replaces each L by L/c and v by the dimensionless preferred velocity v_lab/c. For a=2,v=4/5,L=1, the parallel outbound coordinate time is -1/2 and the incoming coordinate duration is 13/14; the sending clock's total is 3/7>0. A reversed one-way coordinate order does not produce a loop.

Detection would require an actual source/detector coupling, known spectrum and background, latency/noise calibration, and confirmation that rods/clocks retain the assumed matter dynamics. No observed time advance, preferred velocity, experimental sensitivity or exclusion bound is inferred. Existing photon, matter or gravitational-wave limits cannot simply be assigned to this uncoupled scalar.

## Primary-source comparison and reproduction

The targeted live source check consulted:

- Liberati, Sonego and Visser, [gr-qc/0107091v2](https://arxiv.org/abs/gr-qc/0107091v2), revised 14 February 2002, sections 3.1 and 3.2.2: common preferred time and stable causality.
- Bruneton, [gr-qc/0607055v2](https://arxiv.org/abs/gr-qc/0607055v2), revised 2 October 2006, sections II.7, III.3 and V: effective cones, admissible initial surfaces and ultraviolet front limits. Section II.7 explicitly sets fixed non-Lorentz-invariant actions aside; the paper does not validate this toy action.
- Babichev, Mukhanov and Vikman, [0708.0561v1](https://arxiv.org/abs/0708.0561v1), submitted 3 August 2007, sections 4–5 and equations 26–27: a directly comparable homogeneous faster-cone perturbation PDE and its boosted initial-surface problem. Their full k-essence theory differs from this fixed-background action. Section 7's luminal boundary for finite smooth clumps in a regular trivial vacuum has different assumptions; section 6's chronology-protection conjecture is not a proved universal theorem.

[SOURCE_VERIFICATION.json](../research/preferred-frame/SOURCE_VERIFICATION.json) records exact versions, sections, live retrievals and hashes. This is an established-literature comparison, not exhaustive web coverage or a current-discovery search. Article contents are excluded from public records.

From the repository root, with Python>=3.10 and the standard library:

```bash
python3 research/preferred-frame/reproduce.py --out /tmp/ftl-preferred-frame-results.json
python3 research/preferred-frame/independent_check.py --out /tmp/ftl-preferred-frame-independent.json
```

The [registration and amendments](../research/preferred-frame/REGISTRATION.md), [root result](../research/preferred-frame/RESULTS.json), [independent result](../research/preferred-frame/INDEPENDENT_RESULTS.json) and [internal audit](../research/preferred-frame/AUDIT.md) bind the assumptions, strengthened controls, source hashes, outcomes and limitations. The final root run passed 2,000 modes/6,000 derivative components, 2,000 boosts, 21 exact principal matrices, 12 two-resolution energy cells, 2,000 relays, 20 laboratory cells and the boundary/negative controls. The independent route passed 9,206 exact rational checks. These verify conditional methods, not a physical theory's realizability.

Interactions, a dynamical background, quantum/UV consistency, gravitational completion, ordinary relativistic microcausality, production/detection, energy cost and a derived PTA/ringdown link remain unsupplied. Those are requirements for advancing beyond this free illustrative action.
