# Retarded response of the assumed 1+1 dimensional coupled-scalar model

This note derives a response kernel for an explicitly assumed classical scalar system. It is conditional mathematical analysis, using standard retarded Green functions, partial fractions, and Duhamel convolution. It is not a novelty claim, a detected superluminal channel, a coupling to actual ordinary matter, a calibrated detector prediction, or a finite apparatus energy estimate. The numerical verification below is internal and was registered after the algebra was already known.

## Model and units

Set the reference speed to one and use one fixed preferred time. For `a>1`, `m,M>=0` and real `g`, take

\[
\mathcal L=\frac12(\phi_t^2-a^2\phi_x^2-m^2\phi^2+
\chi_t^2-\chi_x^2-M^2\chi^2)-g\phi\chi+J_\chi\chi.
\]

Thus

\[
D_{a,m}\phi+g\chi=0,\qquad D_{1,M}\chi+g\phi=J_\chi,
\quad D_{c,\mu}=\partial_t^2-c^2\partial_x^2+\mu^2.
\]

Only the chi field is sourced and read as an idealized field response. The impulse Green function is a mathematical distribution, not a realizable source or detector specification. Speeds are dimensionless, `m,M` have units of inverse time, and `g` has units of inverse time squared. In these 1+1 dimensional coordinates a free impulse Green function is dimensionless; its two-factor and three-factor convolutions have units of time squared and time to the fourth power. Consequently `g H` and `g^2 K` have the correct Green-function units.

The source-free energy density and flux are

\[
e=\tfrac12(\phi_t^2+a^2\phi_x^2+m^2\phi^2+
\chi_t^2+\chi_x^2+M^2\chi^2+2g\phi\chi),\qquad
j=-a^2\phi_t\phi_x-\chi_t\chi_x,
\]

with `partial_t e + partial_x j = J_chi chi_t`. The assumed stability condition is the positive-semidefinite mass matrix condition

\[
g^2\le m^2M^2.
\]

For spatial wave number k the squared frequencies are

\[
\Omega_\pm^2=\frac{(a^2+1)k^2+m^2+M^2
\pm\sqrt{[(a^2-1)k^2+m^2-M^2]^2+4g^2}}2.
\]

Under the stability condition they are nonnegative. In particular a nonzero stable mixing requires both masses to be positive. The full model with `m=M=0,g!=0` is unstable: at k=0 its squared frequencies are `+|g|` and `-|g|`. Massless convolution coefficients below are useful formal coefficients and majorants; they do not supply a stable nonzero-mixing massless example.

## Retarded convolutions and convergence

The free 1+1 dimensional retarded Green function is

\[
G_{c,\mu}(t,x)=\frac{\theta(t)\theta(ct-|x|)}{2c}
J_0\!\left(\mu\sqrt{t^2-x^2/c^2}\right).
\]

Here star denotes integration over intermediate space and elapsed time. Eliminating the two fields iteratively gives

\[
G_{\phi\chi}=-gH+O(g^3),\quad H=G_{a,m}*G_{1,M},
\]

\[
G_{\chi\chi}=G_{1,M}+g^2K+O(g^4),\quad
K=G_{1,M}*G_{a,m}*G_{1,M}.
\]

More precisely, with `K_n=G_{1,M}^{*(n+1)}*G_{a,m}^{*n}`,

\[
G_{\chi\chi}=\sum_{n\ge0}g^{2n}K_n,\qquad
G_{\phi\chi}=-\sum_{n\ge0}g^{2n+1}
G_{1,M}^{*(n+1)}*G_{a,m}^{*(n+1)}.
\]

For real arguments `|J0|<=1`. At elapsed time s, each free kernel has spatial L1 norm at most s and spatial supremum at most 1/2 for c>=1. Bounding one factor in the supremum norm and the others in L1, and integrating the time simplex, bounds an N-factor convolution by

\[
\frac{t^{2N-2}}{2(2N-2)!}.
\]

In particular `|K_n|<=t^(4n)/[2(4n)!]`. The series converge absolutely and locally uniformly at finite time; the analogous bound also applies to the odd cross series. This convergence holds even without energy stability and must not be confused with it. Every convolution is retarded in the same preferred time and supported within `|x|<=a t`. No finite Fourier cutoff is used to establish this support.

## Explicit massless coefficients

Write `r=|x|`, `delta=a^2-1`, and `q_c=(t-r/c)_+`, with t>=0. The Fourier-Laplace denominators for H0 and K0 are, respectively,

\[
\frac1{(s^2+a^2k^2)(s^2+k^2)},\qquad
\frac1{(s^2+a^2k^2)(s^2+k^2)^2}.
\]

Their decompositions are

\[
\widetilde H_0=\frac1{\delta s^2}
\left[\frac{a^2}{s^2+a^2k^2}-\frac1{s^2+k^2}\right],
\]

\[
\widetilde K_0=\frac{a^4}{\delta^2s^4(s^2+a^2k^2)}
-\frac{a^2}{\delta^2s^4(s^2+k^2)}
-\frac1{\delta s^2(s^2+k^2)^2}.
\]

Use inverse spatial transforms `exp(-s r/c)/(2cs)` and `exp(-s r)(1+s r)/(4s^3)` for the single and double denominators, then invert `exp(-s d)/s^n` as `(t-d)_+^(n-1)/(n-1)!`. This gives

\[
H_0=\frac{a q_a^2-q_1^2}{4\delta},
\]

\[
K_0=\frac{2a^3q_a^4-(3a^2-1)q_1^4-4\delta r q_1^3}
{96\delta^2}.
\]

In the open region between the reference and fast cones, `t<r<a t`, set `Delta=a t-r`. Then

\[
H_0=\frac{\Delta^2}{4a\delta},\qquad
K_0=\frac{\Delta^4}{48a\delta^2}.
\]

The continuous a->1 limits are `theta(t-r)(t^2-r^2)/8` and `theta(t-r)(t^2-r^2)^2/128`. Spatial integration gives `integral H0 dx=t^3/6` and `integral K0 dx=t^5/120`. No fast-front delta occurs in these induced leading kernels; the 1+1 dimensional first-order cross response starts quadratically and the chi response at order g squared starts quartically in Delta.

## Bounded massive triangle integral

At an intercone point let the elapsed slow segment have coordinates `(sigma,y)`. Put `p=sigma+y,q=sigma-y`. Since `r>t`, the fast distance `r-y` is positive throughout the integration domain. The domain is exactly

\[
p,q\ge0,\quad (a-1)p+(a+1)q\le2\Delta.
\]

Use `p=2 Delta u/(a-1), q=2 Delta v/(a+1)` on the unit triangle `u,v>=0,u+v<=1`. The space-time Jacobian is `2 Delta^2/delta`. Define the two nonnegative squared proper intervals

\[
S=\frac{4\Delta^2uv}{\delta},
\]

\[
F=\frac{\Delta}{a}(1-u-v)
\left[2t-\frac{\Delta}{a}
\left(1+\frac{a+1}{a-1}u+\frac{a-1}{a+1}v\right)\right].
\]

Then

\[
H=\frac{\Delta^2}{2a\delta}
\int_{u+v\le1}J_0(m\sqrt F)J_0(M\sqrt S)\,du\,dv.
\]

For the two repeated slow factors,

\[
G_{1,M}*G_{1,M}=-\partial_{M^2}G_{1,M}
=\theta(t-|x|)\frac{\rho J_1(M\rho)}{4M},
\quad \rho^2=t^2-x^2,
\]

with its continuous M=0 limit. Hence, defining `B(z)=2 J1(z)/z` and `B(0)=1`,

\[
K=\frac{\Delta^4}{2a\delta^2}
\int_{u+v\le1}uv J_0(m\sqrt F)B(M\sqrt S)\,du\,dv.
\]

The triangle area is 1/2 and its uv moment is 1/24, recovering H0 and K0. This calculation uses a bounded physical integration domain, rather than a truncated spatial Fourier transform.

For real z, `|J0(z)|<=1`, `|1-J0(z)|<=z^2/4`, `|B(z)|<=1`, and `|1-B(z)|<=z^2/8`. Also `F<=2t Delta/a` and `S<=Delta^2/delta`. Product bounds give

\[
\left|H/H_0-1\right|\le E_H=
\frac{m^2t\Delta}{2a}+\frac{M^2\Delta^2}{4\delta},
\]

\[
\left|K/K_0-1\right|\le E_K=
\frac{m^2t\Delta}{2a}+\frac{M^2\Delta^2}{8\delta}.
\]

Either coefficient is strictly positive when its corresponding error bound is less than one. At fixed t,a,m,M, triangle moments of `1-u-v` give

\[
H/H_0=1-\frac{m^2t\Delta}{6a}+O(\Delta^2),\qquad
K/K_0=1-\frac{m^2t\Delta}{10a}+O(\Delta^2).
\]

These bounds concern the coefficients. A positive order-g-squared coefficient alone is not a proof that the full response is positive.

## Strict positivity and a full-response interval

When `m t<=1` and `M t<=1`, every free segment in a retarded convolution has Bessel argument at most one. The bound `J0(z)>=1-z^2/4` makes all even-mixing terms nonnegative. In the intercone the free chi term is zero, so for g!=0

\[
G_{\chi\chi}\ge
g^2(1-m^2t^2/4)(1-M^2t^2/4)^2K_0>0.
\]

Absolute convergence justifies summing these positive terms. Together with the energy-stability condition, this is a full stable-model mathematical response outside the reference cone. It remains an assumed scalar coupling rather than a physical communication capability.

A sharper remainder estimate follows from the massless repeated-operator kernel

\[
G_{c,0}^{*N}=\theta(t)\theta(ct-r)
\frac{(t^2-r^2/c^2)^{N-1}}
{2c\,4^{N-1}[(N-1)!]^2}.
\]

Collect the n fast and n+1 slow factors and use the same triangle. The moment
`integral u^n v^n(1-u-v)^(n-1) du dv=(n!)^2(n-1)!/(3n+1)!`
and `F<=(2t Delta/a)(1-u-v)` give, for n>=1,

\[
|K_n|\le
\frac{t^{n-1}\Delta^{3n+1}}
{2^n a^n\delta^{n+1}(n-1)!(3n+1)!}.
\]

The n=1 bound equals K0. Define

\[
B_2=\frac{g^4t\Delta^7}{20160a^2\delta^3},\qquad
\eta=\frac{g^2t\Delta^3}{2880a\delta}.
\]

From n=2 onwards successive term bounds have ratio no greater than eta. If eta<1, the absolute sum of all n>=2 terms is at most `B2/(1-eta)`. In the small-time positive regime this provides the rigorous interval

\[
g^2(1-E_K)K_0\le G_{\chi\chi}
\le g^2K_0+\frac{B_2}{1-\eta}
\]

when E_K<1, in addition to the earlier three-factor lower bound. Without small-time positivity, the absolute remainder estimate still applies and gives a lower bound `g^2(1-E_K)K0-B2/(1-eta)`. Consequently, at any fixed positive t and finite parameters, the full chi kernel near the fast front has

\[
G_{\chi\chi}=\frac{g^2\Delta^4}{48a\delta^2}
\left[1+O(m^2t\Delta+M^2\Delta^2+g^2t\Delta^3)\right].
\]

For nonzero g its mathematical support therefore reaches the fast cone. The quartic suppression in this dimension distinguishes support from signal size. The fixed common preferred time prevents these retarded paths from running backwards in that time. None of this identifies actual matter, practical emission, detection, a usable energy budget, or empirical evidence; a physical front in an effective theory would also require a specified ultraviolet completion.

## Verification and correction history

The bounded stdlib verification is specified in REGISTRATION.md and recorded in RESULTS.json. It evaluates 14 massive/front cases covering 12 distinct cells, at Simpson resolutions 64 and 128, plus exact rational transforms, moments and analytic full-response bounds. The computation does not numerically solve or Fourier-invert the complete coupled kernel. The full-response interval follows from the analytic Volterra bounds.

After the initial successful run, the output-path guard was corrected so a public reproduction can write to a fresh writable path. The initial script and results are preserved as INITIAL_REPRODUCE.py and INITIAL_RESULTS.json; their hashes and the original registration hash are recorded in the registration amendment. No mathematical criteria changed. A separate descriptive count error in the original registration, 13 rather than 12 distinct cells, was corrected before the final rerun. The final mathematical payload equals the initial one exactly; only run time and amended source/registration hashes changed. See AUDIT.md for binding and validation details.
