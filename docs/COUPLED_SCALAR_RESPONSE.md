# Local coupling to a preferred fast scalar

Ricardo Maldonado · conditional methods exercise · 2 October 2026 Pacific (3 October UTC)

The [free preferred-time scalar exercise](PREFERRED_FRAME_SCALAR.md) supplied a propagation equation but no source/detector interaction. This note adds a **specified local scalar analogue**. A stable parameter region permits an idealized source acting on the slower scalar to produce a nonzero response of that same scalar outside its uncoupled light cone. The response is suppressed near the faster cone. These are familiar linear-field methods applied to an assumed model, not a novel theory, physical FTL demonstration, measured signal, completed hidden-sector model or external peer review.

## Assumed model and energy

Use ordinary c=1 and the previous fixed global preferred time. Let phi have principal speed a>1 and chi principal speed1. Assume

$$
\mathcal L=\tfrac12[\dot\phi^2-a^2|\nabla\phi|^2-m^2\phi^2
+\dot\chi^2-|\nabla\chi|^2-M^2\chi^2]-g\phi\chi+J\chi,
\qquad m,M\ge0.
$$

The real mixing g has units of mass squared. J is prescribed classical forcing and chi is the mathematical response; neither is an actual photon/matter emitter or detector. The equations are

$$
(\partial_t^2-a^2\Delta+m^2)\phi+g\chi=0,
\qquad(\partial_t^2-\Delta+M^2)\chi+g\phi=J.
$$

The preferred internal-field energy density and flux are

$$
e=\tfrac12[\dot\phi^2+\dot\chi^2+a^2|\nabla\phi|^2+|\nabla\chi|^2
+m^2\phi^2+M^2\chi^2]+g\phi\chi,
\quad\mathbf S=-a^2\dot\phi\nabla\phi-\dot\chi\nabla\chi,
\quad\partial_t e+\nabla\cdot\mathbf S=J\dot\chi.
$$

For fields on unrestricted R^d, or periodic boxes admitting k=0, the source-free Hamiltonian is bounded below precisely when the mass matrix is positive semidefinite. A fixed bounded Dirichlet domain can instead have a stabilizing gradient gap; that restricted-domain exception is outside the k=0 criterion used here. The externally forced density written above has canonical Hamiltonian `e-J*chi`; the work equation uses the internal field energy e.

$$\boxed{g^2\le m^2M^2.}$$

For m>0, the potential is `(m*phi+g*chi/m)^2/2+(M^2-g^2/m^2)*chi^2/2`. At equality one uniform mass eigenmode is zero where such modes are allowed. If either mass is zero, stable mixing requires g=0. In particular, **two massless scalars with nonzero mass mixing are unstable**; massless kernels below describe perturbative coefficients, not a stable massless coupled system. Positive free energy does not supply the energy, bandwidth or backreaction of a realizable source apparatus.

## Modes and the two characteristic cones

Set A=a^2|k|^2+m^2, B=|k|^2+M^2 and D=sqrt((A-B)^2+4g^2). Then

$$\omega_\pm^2=\tfrac12(A+B\pm D).$$

Under the stability condition all preferred frequencies are real, and the matrix after subtracting |k|^2 times the identity remains positive semidefinite, so omega_±^2>=|k|^2. A g^2>m^2M^2 low-k tachyon is a potential instability; a wrong kinetic sign is a different ghost pathology.

The fixed background breaks physical boost symmetry, although coordinates and the background can be transformed together. Preferred-positive modes have positive ordinary-boost frequencies because omega>=|k|. A boosted t'=constant plane is jointly spacelike for both principal cones only when |v|*a<1. Apparent arbitrary-data growth on a non-Cauchy plane is a wrong-slicing obstruction, not an additional preferred-frame mass or kinetic instability; the [free-field discussion](PREFERRED_FRAME_SCALAR.md#boosted-initial-surfaces-and-apparent-negative-coefficients) distinguishes the corresponding charges and fluxes.

Local mass mixing does not alter the principal symbol. The system retains both speed1 and speeda characteristic cones, with outer causal support bounded by the speeda cone. Its retarded solution uses one increasing preferred time. Every term of the causal Volterra expansion stays within that outer cone. Common-time retarded channels and ordinary future-directed relays cannot form a finite loop under these assumptions. This does not extend to unspecified backgrounds, quantum/gravitational completions or other channels that violate the common ordering.

In the full continuum model the chi-to-chi response has the outer front a when g!=0, even though chi's isolated principal equation has speed1. At high |k|, the fast normal branch has chi spectral weight

$$w_{\chi,+}\sim\frac{g^2}{(a^2-1)^2|k|^4}.$$

Thus the fast contribution is suppressed, rather than absent. A band-limited Fourier plot cannot establish a sharp causal front. An effective low-energy action alone also cannot determine the physical ultraviolet front. Momentum-dependent spectral projections are spatially nonlocal; exterior tails of a separately projected normal branch do not contradict the support of the complete local coupled response.

## Retarded response in one spatial dimension

Normalize the free kernels by `(partial_t^2-c^2*partial_x^2+mu^2)G=delta(t)delta(x)`. In 1+1 dimensions,

$$G_{c,\mu}^R(t,x)=\frac{\theta(t)\theta(ct-|x|)}{2c}
J_0\!\left(\mu\sqrt{t^2-x^2/c^2}\right).$$

Write star for spacetime convolution, G_f=G_{a,m}^R and G_s=G_{1,M}^R. Retarded elimination gives

$$G_{\phi\chi}=-gG_f*G_s+O(g^3),\qquad
G_{\chi\chi}=\sum_{n=0}^\infty g^{2n}G_s^{*(n+1)}*G_f^{*n}.$$

The cross-field response starts linearly in g; sourcing and reading chi requires two conversions and starts quadratically outside the slower cone. No source-energy normalization or finite detector response is inferred from a delta-source kernel.

For the massless coefficient functions, let r=|x|, q_c=max(t-r/c,0), d=a^2-1, with kernels zero for t<=0. Exact convolution gives

$$H=G_{a,0}^R*G_{1,0}^R=\frac{a q_a^2-q_1^2}{4d},$$

$$K=G_{1,0}^R*G_{a,0}^R*G_{1,0}^R
=\frac{2a^3q_a^4-(3a^2-1)q_1^4-4drq_1^3}{96d^2}.$$

They are normalized by `integral(H dx)=t^3/6` and `integral(K dx)=t^5/120`. In the intercone region t<r<a*t,

$$\boxed{H=\frac{(at-r)^2}{4a(a^2-1)},\qquad
K=\frac{(at-r)^4}{48a(a^2-1)^2}.}$$

The ordinary free response vanishes there. These coefficients vanish at and beyond the fast front. At a=1 use the separate coincident-speed kernels `H=(t^2-r^2)/8` and `K=(t^2-r^2)^2/128` inside r<t, rather than substituting into singular distinct-speed expressions.

## A stable massive response, rather than an unstable massless example

For small positive t with m*t<=1 and M*t<=1, each constituent massive kernel has `1-mu^2*t^2/4<=J0<=1` on its causal support. All even-mixing terms in the chi-to-chi Volterra series are then nonnegative. The term with 2n+1 factors is bounded in absolute value by `t^(4n)/(2*(4n)!)` for a>=1, which ensures local absolute convergence after multiplying by g^(2n). In t<r<a*t the free slower term is zero, and

$$G_{\chi\chi}(t,r)\ge
g^2\left(1-\frac{m^2t^2}{4}\right)
\left(1-\frac{M^2t^2}{4}\right)^2
\frac{(at-r)^4}{48a(a^2-1)^2}>0$$

for g!=0. This bound is specifically **1+1 dimensional**. Choose a=2,m=M=1,g=1/2 as a stable example. It proves a nonzero **mathematical** same-field early response without identifying it with a physical signaling device. It does not estimate detectable signal-to-noise or engineering feasibility.

Mass corrections to the leading fast-front coefficient vanish with at-r; further mixing terms are more suppressed there. The kernel review separately derives a bounded triangle quadrature for the massive order-g^2 term. Its numerical checks are checks of that coefficient, not direct numerical evaluation of the complete Volterra sum.

## Spatial dimension and effects on the slower scalar

For a rotationally invariant 3+1 kernel at r>0, componentwise dimensional reduction gives `G^(3)=-(2*pi*r)^-1*partial_r G^(1)`. Consequently the massless intercone coefficients are

$$H_3=\frac{at-r}{4\pi a(a^2-1)r},\qquad
K_3=\frac{(at-r)^3}{24\pi a(a^2-1)^2r}.$$

The induced same-field onset is cubic in this radial coefficient, not the one-dimensional quartic. The relation must include distributions at ordinary fronts; r=0/contact terms are outside these explicit formulas. The 3+1 free massive kernel has a delta front and a negative interior Bessel tail, so differentiating the 1+1 inequality does not establish a 3+1 positive-response bound. These are mathematical point-source coefficients, not experimental predictions.

The coupling also changes the slower field's dispersion. For m^2>M^2 the slope of its lower branch at k^2=0 is

$$\left.\frac{d\omega_-^2}{d|k|^2}\right|_0
=\frac{a^2+1}{2}-\frac{(a^2-1)(m^2-M^2)}{2\sqrt{(m^2-M^2)^2+4g^2}}>1\quad(g\ne0).$$

For a massive mode this coefficient is **not** its k=0 group speed or physical front. Retarded elimination gives the exact chi operator `D_s-g^2*(D_f^R)^-1`. If the phi mass is heavy and `|omega^2-a^2*k^2|<<m^2`, its derivative expansion has `Z_t=1+g^2/m^4`, `Z_x=1+a^2*g^2/m^4`, and unnormalized mass coefficient `M^2-g^2/m^2`. After canonical normalization, the truncated mass squared is `(M^2-g^2/m^2)/Z_t` and the gradient coefficient is `Z_x/Z_t`. Higher derivative terms remain; no UV-front prediction follows from truncation. It is therefore inconsistent to identify chi with otherwise unchanged ordinary matter and carry over ideal rods/clocks/timing assumptions without deriving that matter model.

## Verification and remaining work

The [registered reproduction and internal reviews](../research/coupled-scalars/README.md) bind exact algebra, modes, support coefficients, normalization, stable massive bounds, and pathological controls. They are synthetic internal assistance, not observations or external replication. The initial independent control labels and a kernel output-path restriction were corrected transparently; the final result binds the strengthened sources.

The preceding [finite delay](FINITE_RESPONSE_DELAY.md), [unequal sender-speed](ASYMMETRIC_SIGNAL_SPEEDS.md) and [free preferred-time](PREFERRED_FRAME_SCALAR.md) assumptions remain distinct. Common-time coupling does not revive the original sender-frame loop diagnostic, nor supply a calibrated PTA/ringdown connection. See the prior [primary-source comparison](PREFERRED_FRAME_SCALAR.md#primary-source-comparison-and-reproduction) for its exact-version preferred-causality and front-scope citations; this round adds a derivation, not a current-discovery literature claim.

Open requirements include physical field identification, controllable finite-energy source/detector dynamics and calibration, constraints on transferred Lorentz violation, quantum/UV and gravitational consistency, source backreaction, realistic noisy predictions, and recovery of missing private v6.3 code/proofs and authenticated strict observational inputs. This exercise supplies an assumed coupling and a rigorous response calculation; it does not settle those requirements.
