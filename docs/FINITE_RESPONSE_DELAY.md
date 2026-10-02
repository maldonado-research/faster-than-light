# Finite proper reply delay in the FTL antitelephone model

Ricardo Maldonado · subsequent methods exercise · 2026-10-01 Pacific (2026-10-02 UTC)

The archived construction assumes Bob replies immediately. This note changes that assumption explicitly and derives its consequences. It is conditional kinematics, not a demonstrated communication channel, observational result, external peer review, or claim of mathematical novelty. It does not alter the September 8 snapshot or archived equations. The baseline and its assumptions are in [CORE_MATHEMATICS.md](CORE_MATHEMATICS.md) and [RESEARCH_ACCOUNT_AND_EQUATIONS.md](RESEARCH_ACCOUNT_AND_EQUATIONS.md).

## Model and event calculation

Alice and Bob coincide at time zero. In Alice's frame Bob moves at speed beta*c, with 0<=beta<1. Alice sends at proper time T>0. Each sender can hypothetically transmit at speed W*c, W>1, in their own frame. Motion is collinear and inertial in Minkowski spacetime. Bob waits a **proper** interval tau>=0 after reception before replying.

Define

$$
s=\sqrt{1-\beta^2},\qquad d=\frac{\tau}{T},\qquad
A=\frac{Ws}{W-\beta}.
$$

The outgoing worldline intersection is t1=WT/(W-beta), x1=beta*c*t1. Lorentz transformation gives Bob's reception time t1'=A*T. He departs at t_r'=A*T+tau from x'=0. In Bob's frame the return ray x'=-W*c*(t'-t_r') intersects Alice's worldline x_A'=-beta*c*t' at

$$
t_2'=\frac{W(AT+\tau)}{W-\beta}.
$$

Alice's return event has x=0, so t2=s*t2'. Therefore

$$
\boxed{\frac{t_2}{T}=A^2+dA
=\frac{W^2(1-\beta^2)}{(W-\beta)^2}
+d\frac{W\sqrt{1-\beta^2}}{W-\beta}.}
$$

The Alice-frame duration of Bob's wait is gamma*tau. Adding tau directly to Alice's reception time would confuse proper and coordinate time.

## Delay allowance and exact threshold

For fixed W and beta, return before sending requires

$$
d<d_{\max}=A^{-1}-A.
$$

A nonnegative allowed delay exists only inside the immediate-reply region beta>2W/(1+W^2). Equality d=d_max gives t2=T. For example W=2 and beta=12/13 give A=5/7 and d_max=24/35. At d=1/2, t2/T=85/98<1.

For fixed finite d>=0, put

$$
h=\frac{2}{\sqrt{d^2+4}+d}.
$$

Since A>0, A^2+dA<1 exactly when A<h. Solving gives

$$
\boxed{\beta>\beta_d(W)
=\frac{W\left[h^2+\sqrt{W^2(1-h^2)+h^2}\right]}{W^2+h^2}.}
$$

This recovers 2W/(1+W^2) at d=0. Every finite d leaves beta_d<1, so a finite proper delay does not remove the modeled region for all observer speeds. Close to beta=1, evaluate the equivalent gap to avoid cancellation:

$$
1-\beta_d=
\frac{h^2(W-1)^2}
{W^2-h^2(W-1)+W\sqrt{W^2(1-h^2)+h^2}}.
$$

The return ratio increases with delay and decreases with W for beta>0. It is **not globally decreasing in beta**:

$$
\partial_\beta A=
\frac{W(1-W\beta)}{\sqrt{1-\beta^2}(W-\beta)^2}.
$$

It initially increases, reaches a maximum at beta=1/W, and then decreases. This does not obstruct the unique positive-beta threshold on the decreasing branch.

## Observer-speed caps and timing assumptions

If beta is restricted to [0,beta_max] with beta_max<1 and beta_max>2W/(1+W^2), define A_cap=A(beta_max). No return before sending occurs anywhere in that interval exactly when

$$
d\ge A_{\rm cap}^{-1}-A_{\rm cap}.
$$

Equality gives t2=T at the cap. If beta_max is already outside the immediate-reply region, a delay is unnecessary for this limited construction.

The conclusion holds at fixed d=tau/T. A fixed finite absolute tau does not uniformly suppress the modeled region when T can grow: d approaches zero. No result here establishes that a delay removes every causal loop in a dynamical theory.

## The population ceiling changes

Define the delayed positive strength

$$
P_{d,+}=\max(1-A^2-dA,0),\qquad
J_d(W)=\mathbb E[P_{d,+}^2].
$$

At fixed d and a W-independent beta population,

$$
\boxed{J_{\infty}(d)=
\mathbb E\left[(\beta^2-d\sqrt{1-\beta^2})_+^2\right].}
$$

The archived E[beta^4] is recovered at d=0; it cannot be carried unchanged into this delayed model. If delays depend on beta, W or T, their population distribution must be specified separately.

For **uniform beta** and fixed d, the infinite-W boundary is beta0=sqrt(1-h^2). The change beta=cos(u) yields the smooth integral

$$
J_\infty(d)=\int_0^{\arcsin h}
(\cos^2u-d\sin u)^2\sin u\,du.
$$

The sin(u) factor preserves the uniform-beta measure. A separate positive representation, using sqrt(1-beta^2)=h*z and d*h=1-h^2, is

$$
J_\infty(d)=h^2\int_0^1
\frac{z(1-z)^2(1+h^2z)^2}{\sqrt{1-h^2z^2}}\,dz.
$$

It follows that the ceiling is positive for every finite d, strictly decreases with d, and satisfies

$$
J_\infty(0)=\frac15,\qquad
\lim_{d\to\infty}d^2J_\infty(d)=
\int_0^1z(1-z)^2\,dz=\frac1{12}.
$$

| Proper delay ratio d=tau/T | Uniform-beta ceiling J_infinity(d) |
|---:|---:|
| 0 | 0.200000000000000 |
| 0.1 | 0.165751266498482 |
| 0.5 | 0.086389378194352 |
| 1 | 0.044529385232735 |
| 2 | 0.016655437895056 |
| 10 | 0.000824309012249 |
| 100 | 0.000008332416812 |
| 1000 | 0.000000083333242 |

If one chooses to retain the author's proposed x=20*DeltaNeff*J diagnostic in this extended model, the finite-crossing criterion must use J_infinity(d), with a strict inequality J_bound<J_infinity(d). Neither this replacement nor the original diagnostic derives an FTL interaction, energy requirement or calibrated discovery threshold.

## Reproduction and internal audit

The [registered protocol](../research/response-delay/REGISTRATION.md) discloses prior algebraic and agent knowledge. Run the standard-library implementation from the repository root:

```bash
python3 research/response-delay/reproduce.py --out /tmp/ftl-response-delay-results.json
```

Python 3.12.14 was used; Python>=3.10 is required. No package, service, secret or observation is needed. The command preserves failure through a nonzero exit status. The [saved result](../research/response-delay/RESULTS.json) records parameters, check counts, residuals and source hash; tolerances are specified in the registration and implementation.

The current run passed 5,000 direct event comparisons, 72 exact rational event comparisons, 30 high-precision thresholds and 60 boundary-side checks. Eight population ceilings were checked at two resolutions, together with 30 finite-W severity cells. Maximum event relative difference was 2.60e-14, threshold residual was 1e-69, and relative ceiling resolution change was 7.22e-15.

A separate internal agent derived the rescaled population integral and independently checked thresholds, caps and ceilings. This is internal AI-assisted mathematical checking, not external peer review. See [AUDIT.md](../research/response-delay/AUDIT.md). No physical FTL signal, production process, stability result or empirical constraint is established.
