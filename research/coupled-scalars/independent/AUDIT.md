# Independent internal audit: local coupled-scalar analogue

3 October 2026 UTC / 2 October Pacific. This is an explicitly assumed
quadratic classical model and internal mathematical verification. It is not
a physical FTL channel, Ricardo's established hidden-sector theory, a detector
proposal, novelty, observation, external peer review or empirical replication.
No root implementation was read before the registered independent run.

## Model, Hamiltonian and stability

Write mu=m²>=0, nu=M²>=0, a>1 and real g. The stipulated equations are

    (partial_t²-a² Delta+mu)phi+g chi=0,
    (partial_t²-Delta+nu)chi+g phi=J_chi.

The preferred momenta are pi_phi=phi_t and pi_chi=chi_t. The free energy is

    e=1/2[pi_phi²+pi_chi²+a²|grad phi|²+|grad chi|²]
        +1/2[mu phi²+2g phi chi+nu chi²],
    S=-a² pi_phi grad phi-pi_chi grad chi,
    partial_t e+div S=J_chi*pi_chi.

Here e is the free/internal field energy. The canonical Hamiltonian density
of the Lagrangian with prescribed forcing is e-J_chi*chi; arbitrary forcing
is not covered by the source-free energy-positivity claim.

The potential matrix V=[[mu,g],[g,nu]] is positive semidefinite exactly
when g²<=mu*nu. With the already assumed positive kinetic and gradient
terms this is the all-momentum stability and Hamiltonian lower-bound criterion
on flat unbounded space (or when a periodic zero mode is admitted). Necessity
can also be seen using wide localized fields in a negative mass direction:
their potential dominates the gradients and scaling their amplitude makes H
unbounded below. A finite box with a strictly positive spatial eigenvalue gap
would require its own spectral boundary analysis.

For mu>0 the potential square completion is

    2V_potential=(m phi+g chi/m)²+(nu-g²/mu)chi².

At g²=mu*nu a k=0 null direction exists. A constant zero-energy field is
admitted only under appropriate boundaries; decaying whole-space states retain
a gapless small-k branch. If either mass is zero, stable mixing requires g=0.
The signs +/-g have identical spectra and same-chi response, while cross-field
response changes sign. Stable positive/negative saturation and zero-mass
boundaries are explicitly included in the independent checks.

At q=|k|² define A=a²q+mu, B=q+nu and D=sqrt[(A-B)²+4g²]. The modes are

    lambda_+/-=omega_+/-²=(A+B+/-D)/2.

The eigenvectors are momentum dependent. If g²>mu*nu, lambda_-(0)<0:
this is a low-momentum tachyonic instability rather than a wrong-sign kinetic
ghost. Its unstable band is q<q_c, with

    q_c=[-(mu+a²nu)+sqrt((mu-a²nu)²+4a²g²)]/(2a²).

A reversed preferred kinetic sign is a separate ghost control; a reversed
gradient sign is a separate high-k instability. None of them is repaired by
calling an inadmissible boosted slice a new physical state.

For stable V, K(q)-qI=V+diag((a²-1)q,0) is positive semidefinite. Hence
omega_+/->=|k| and ordinary boosted positive-frequency mode energies
gamma(omega-v k_x) stay positive. More generally, physical momentum is
P_x=-integral(pi_phi phi_x+pi_chi chi_x), and the preferred-slice charge is
Q_(partial_t')=gamma(H-vP_x). Pointwise its unscaled integrand equals

    1/2[(pi_phi+v phi_x)²+(pi_chi+v chi_x)²
      +(a²-v²)phi_x²+(1-v²)chi_x²
      +a²|grad_perp phi|²+|grad_perp chi|²]+V_potential >=0

for |v|<1. The universal bound is H>=|P_x|, governed by the slow sector,
rather than the previous single-fast-field bound H>=a|P_x|. Prescribed forcing
can supply or remove energy; these conserved-charge statements concern J=0.

## Cones, observers and ultraviolet scope

Zero-order mass/mixing leaves two principal characteristic factors,

    (tau²-a²|xi|²)(tau²-|xi|²).

There are two characteristic cones and one outer domain-of-dependence cone
of speed a. For the stable energy, |S dot n|<=a e. The moving-boundary energy
argument therefore excludes response outside the fast cone for retarded data
and sources. Duhamel support below proves this directly as well, even before
using a positivity criterion. A tachyonic system can remain causal in this
mathematical sense while failing the stability benchmark.

The boosted fast kinetic coefficient is gamma²(1-a²v²), while the ordinary
chi coefficient is 1. Common spacelike causal Cauchy surfaces require
|v|a<1. Equality is characteristic; above it the fast surface is non-Cauchy.
At large transverse momentum, a negative fast coefficient gives unbounded
imaginary-frequency growth for arbitrary data on that surface; finite algebraic
mixing cannot alter the principal behavior. The earlier 1+1 caveat still
applies. A negative formal t'-slice kinetic term is not by itself a ghost of
this preferred-time stable theory; translation charges are evaluated on valid
preferred Cauchy surfaces.

At high momentum, Delta=A-B grows as (a²-1)q. Then

    lambda_+=a²q+mu+O(q^-1),
    lambda_-=q+nu+O(q^-1),
    fast chi weight=(1-Delta/D)/2
        =2g²/[D(D+Delta)] ~ g²/[(a²-1)²|k|⁴].

This is mathematical suppression of fast-sector modal weight, not a detector
sensitivity or cutoff. Complete-continuum front statements assume the full
operator remains valid at arbitrarily high frequency. An unspecified EFT
ultraviolet completion does not fix a physical front.

Every actual preferred-retarded leg and ordinary future relay increases the
same t. A finite nontrivial chain cannot close or return to an earlier event
on its sender's worldline. This assumes all channels and source prescriptions
respect the same preferred temporal ordering; unspecified future-dependent
sources, nonlocal interactions or variable backgrounds require new analysis.

## Exact retarded matrix and channel order

For each real k the stable retarded Fourier kernel is

    R(t,k)=theta(t) sin[t sqrt(K(k))]/sqrt(K(k)),

with a zero eigenvalue interpreted by sin(t sqrt(lambda))/sqrt(lambda)->t.
Spectral projectors provide an exact matrix representation, verified independently
against the matrix power series. The full space-time Green matrix is its
distributional inverse Fourier transform.

Let F=G_a,m and C=G_1,M denote the two **free** preferred-retarded kernels,
and let star include space and time convolution. Eliminating phi causally gives

    G_chichi=C+sum_(n>=1) g^(2n) C^*(n+1)*F^*n,
    G_phichi=-sum_(n>=0) g^(2n+1) F^*(n+1)*C^*(n+1).

The first cross-field response is -g F*C, but the first extra chi-source,
chi-measured response is g² C*F*C. A chi-only measurement therefore does not
get a first-order-in-g fast signal. Changing g's sign changes the cross channel
and leaves the same-chi channel unchanged.

Each convolution term is retarded, and its support is a Minkowski sum of
cones whose speeds are <=a. Thus each term vanishes outside |x|<=a t and
the convergent full series does too. Momentum-dependent normal-branch
projectors are nonlocal spatial operators. Their separated kernels need not
have independent compact cones; exterior tails in such a decomposition cancel
in the complete local-PDE response. They cannot be treated as autonomous
signaling channels. A finite Fourier cutoff produces tails and is not a proof
of a sharp causal front.

## Closed massless coefficients, with the instability warning

In 1+1 dimensions G_c,0=theta(t)theta(c t-|x|)/(2c). Set r=|x|, d=a²-1,
H=G_a,0*G_1,0 and K=G_1,0*G_a,0*G_1,0. Both vanish for t<=0 or r>=a t.
For 0<=r<=t,

    H=(a t²-r²)/[4a(a+1)],
    K=[a(2a+1)t⁴-6a t²r²+(a+2)r⁴]/[96a(a+1)²].

For t<r<a t, with delta=a t-r,

    H=delta²/[4a d],
    K=delta⁴/[48a d²].

An independent spatial route uses H''=(G_a,0-G_1,0)/d and
K''=(H-G_1,0*G_1,0)/d, with evenness and zero exterior value/slope.
The exact normalizations are integral(H dx)=t³/6 and integral(K dx)=t⁵/120.
At a=1 they become (t²-r²)/8 and (t²-r²)²/128 inside the common cone.

The full m=M=0,g!=0 model is unstable. These massless expressions are
perturbative coefficients, not a stable massless propagation mechanism.

## Nonzero response in the full stable massive system

In 1+1 dimensions,

    G_c,m=theta(t)theta(c t-r) J0(m sqrt(t²-r²/c²))/(2c).

For 0<=z<=1 the alternating Bessel series gives 1-z²/4<=J0(z)<=1.
If max(m,M)T<=1, every subinterval in a retarded convolution with total
t<=T satisfies this bound. Hence all free kernels are nonnegative on that
domain, and

    G_a,m >= b_phi G_a,0,    b_phi=1-m²T²/4>0,
    G_1,M >= b_chi G_1,0,    b_chi=1-M²T²/4>0.

All same-chi Volterra terms carry even powers of real g and are nonnegative.
The zeroth C term vanishes in t<r<a t, but the first mixed term is strictly
positive there. For g!=0,

    exact G_chichi(t,r) >= g² b_chi² b_phi K(t,r) >0.

This certifies a nonzero **exact full stable response**, rather than merely
an unstable massless Born coefficient. Convergence follows from choosing one
fast factor with sup norm 1/(2a) and bounding the spatial L1 norm of each other
factor by its elapsed time. The nth mixed term is bounded by

    |g|^(2n) t^(4n)/[2a(4n)!], n>=1.

The entire majorant converges for every finite t. For the positivity regime
its positive sum is also an upper bound; RESULTS records exact positive lower
bounds and a six-term majorant with a rational geometric tail bound.
The global bound |J0(z)|<=1 for real z follows from its angular cosine integral;
it supplies absolute convergence beyond the restricted small-time positivity
regime. Positivity itself is claimed only in the stated z<=1 regime.

The registered 54 examples use a in {3/2,2,3}, m=M=1, g=+/-1/2, three
small times and three strict intercone positions. They satisfy the positive-
definite mass criterion. All lower bounds are strictly positive. The script
computes certificate identities and bounds; it does not report a numerical
measurement of the full Green function or an experimental signal amplitude.
At longer times massive Bessel kernels can change sign; no claim that the
response is positive or nonzero at every intercone point is made.

## Front power and formal 3D dimension raising

For fixed finite parameters and t>0, approaching the fast front from its
interior with delta=a t-r down to zero gives in 1+1 dimensions

    G_phichi=-g delta²/[4a d]+O(delta³),
    G_chichi=g² delta⁴/[48a d²]+O(delta⁵).

The slow free term is absent there. The leading coefficients are independent
of diagonal masses: J0 tends to one on a fast characteristic. A fast mass
insertion adds at least one extra delta power; an additional slow kernel adds
two; a higher mixing pair adds at least three. The locally convergent series
therefore cannot cancel the displayed leading terms. The exact support reaches
the fast front for nonzero mixing, although the same-chi onset is soft and the
value at the boundary itself vanishes.

For a rotationally invariant component of this same constant-coefficient PDE,
the distributional dimension-raising identity at r>0 is

    G^(3)(t,r)=-(2*pi*r)^(-1) partial_r G^(1)(t,r).

It follows directly from the 1D and 3D radial Fourier inversions. For the
massless convolution coefficients in the intercone,

    H3=delta/[4*pi*a*d*r],
    K3=delta³/[24*pi*a*d²*r].

For the full massive coupled system, the corresponding leading cross/same-chi
orders are -g H3+O(delta²) and g² K3+O(delta⁴). Thus a quartic **1D**
same-chi onset becomes cubic in this formal **3D** radial kernel. The r>0
formula is not evaluated at the origin, and causal/distributional boundaries
are retained. Exact 3D radial integrals reproduce t³/6 and t⁵/120. These are
mathematical kernel identities, not predictions for a physical source/detector
apparatus in three dimensions.

## Mixing does not preserve an unchanged ordinary matter sector

Chi's uncoupled principal speed is one; with mixing, its response has access
to the faster characteristic. It is an ordinary-principal-speed scalar analogue,
not derived real photon/matter dynamics or validated unchanged rods/clocks.
For mu>nu and g!=0, the chi-like low-k branch has

    d lambda_-/d q at 0 = (a²+1)/2
      -(a²-1)(mu-nu)/(2 sqrt[(mu-nu)²+4g²]) >1.

This is a dispersion gradient, not a front; the massive group velocity still
approaches zero at k=0. At nonzero saturation with both masses positive, the
gapless branch has speed squared (a²nu+mu)/(mu+nu), between one and a².

Causal heavy-phi elimination gives D_chi-g² D_phi,ret^(-1). For
p=omega²-a²k² with |p|<<mu, the local derivative expansion has

    Zt=1+g²/mu², Zx=1+a²g²/mu², mass coefficient C0=nu-g²/mu,
    exact-minus-truncated symbol=-g²p²/[mu²(mu-p)].

After canonical normalization, the truncated effective mass squared is
C0/Zt and its gradient-speed coefficient is Zx/Zt. C0 alone is the
unnormalized mass coefficient, not the normalized mass squared.
The remainder was tested exactly in 60 rational cells. Extrapolating this
truncation to ultraviolet frequency would violate its stated regime. Its
coefficient ratio cannot be assigned to actual clock matter or used as a
physical-front measurement.

## Reproduction and limits

Run:

```bash
python3 /workspace/ftl-rounds/20261002-round-05/independent/check.py \
  --out /tmp/ftl-coupled-independent.json
```

The strengthened final run passed 10,415 checks using Python 3.12.14 and the
standard library, with 54 positive full-response certificates. Exact fraction
equality and independent matrix-series/spectral comparisons were used.
RESULTS.json saves the script and registration hashes, parameter cases,
maximum residuals, actual injected candidates, controls and empty failure list.
Mathematical stability,
support and a nonzero scalar-field response are established conditionally;
production energy, finite-energy pulses, actual detector observables, ordinary
matter consistency, UV completion, interactions beyond quadratic mixing,
gravity, empirical calibration and an observational PTA/ringdown bridge are
not supplied. No private raw material, third-party article contents or claimed
physical observations were used or copied.

## Correction history and strengthened controls

The first run reported 10,414 passes. Skeptical review identified that its four
purported mutation controls merely compared correct formulas with wrong values
or an expected unstable condition; they did not replace a candidate and feed
it through a checker. Those four initial pass flags are not certified mutation
rejection. The initial negative-gradient and negative-kinetic controls also
used literal negative expectations rather than shared parameterized helpers.
The initial mathematical derivations and other checks are retained separately;
these control gaps are disclosed rather than hidden by a new result file.

INITIAL_CHECK.py, INITIAL_RESULTS.json, INITIAL_AUDIT.md, INITIAL_MANIFEST.json
and INITIAL_REGISTRATION.md preserve the original artifacts. The independent
registration now contains a dated post-run amendment explaining the weakness,
the replacement criteria and the anticipated count change. This correction
is not blind preregistration. No root implementation was read for either run.

The final code injects four actual candidate-function defects in memory into
the same 63-constraint property checker used by the healthy candidate:

- The healthy candidate passes all 63 constraints.
- A mass predicate accepting massless nonzero mixing fails one constraint.
- An odd-in-g same-chi response fails 12 constraints.
- A nonzero exterior response fails 12 constraints.
- Copying the 1D same-chi kernel into the radial 3D response fails 24 constraints.

The gate checks actual returned candidate classifications/responses against
low-k spectral behavior, parity, scaling, future support, dimension raising
and radial front power. Each injected replacement and its rejection reasons
is saved in RESULTS. Inner gate evaluations are reported separately; they are
not added to the headline total. Replacing the same four initial outcomes and
adding one healthy-candidate outcome raises the total from 10,414 to 10,415.

The final negative-gradient control evaluates the shared dispersion helper
at reversed and healthy gradient parameters and compares actual least-mode
eigenvalues. It also evaluates an energy state through the shared energy
helper. The negative-kinetic control evaluates the same velocity state through
that energy helper with kinetic signs -1 and +1. The same helpers serve the
healthy original cells. The two old outcome slots are strengthened rather
than counted twice.

The final script requires --out, refuses an existing destination or symlink
before executing checks, and creates results exclusively. After the single
strengthened scientific run, a separate invocation from /tmp attempted the
existing RESULTS.json. It exited 2 with the explicit refusal and unchanged
result SHA-256; OVERWRITE_GUARD.json records this CLI check. This invocation
did not execute the scientific suite and is not included in its count.
Registration and data paths derive from __file__, so the root/independent
layout can be retained under a central registry without depending on cwd.

A separate code-only internal reviewer checked the strengthened candidate
gate, shared parameterized helpers, output guard and correction disclosure
and found no material concern. That reviewer ran no code and added no pass
count. The gate demonstrates sensitivity to the four named defects; it is
not a certificate of exhaustive candidate correctness.
