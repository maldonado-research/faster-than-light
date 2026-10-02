# Unequal signal speeds with a proper reply delay

Ricardo Maldonado · subsequent methods exercise · 2026-10-01 Pacific (2026-10-02 UTC)

The equal-speed antitelephone assumption can be relaxed without supplying a physical FTL mechanism. This note derives the unequal-speed result, identifies an exact loop/no-loop counterexample when the speeds are swapped, and states the assumptions under which a preferred time would exclude the loop. These are conditional mathematical results with internal checking, not a novelty claim, physical signal observation or external peer review. The historical research snapshot remains preserved.

## Events and assumptions

Alice and Bob coincide at time zero and move inertially in Minkowski spacetime. Alice sends at T>0; Bob recedes at beta*c, 0<=beta<1. The hypothetical controllable outbound speed is U*c in Alice's frame and the return speed is V*c in Bob's frame. Initially U,V>1. Bob waits proper time tau>=0 after reception; d=tau/T. Directions are collinear.

Set s=sqrt(1-beta^2), A=U*s/(U-beta), B=V*s/(V-beta). With c=1, the outgoing ray meets Bob at t1=U*T/(U-beta), x1=beta*t1. Bob's reception time is t1'=A*T. The return ray, emitted at t_r'=A*T+tau, intersects Alice's worldline x_A'=-beta*t' at t2'=V*t_r'/(V-beta). The reception event has x=0 in Alice's frame, hence t2=s*t2'. Therefore

$$
\boxed{R\equiv\frac{t_2}{T}=B(A+d)
=\frac{UV(1-\beta^2)}{(U-\beta)(V-\beta)}
+\frac{dV\sqrt{1-\beta^2}}{V-\beta}.}
$$

U=V=W recovers the [preceding proper-delay result](FINITE_RESPONSE_DELAY.md). The delay is Bob's proper interval; adding tau directly in Alice's frame would be incorrect.

## Immediate and delayed thresholds

For d=0,

$$
1-R_0=\frac{\beta[(UV+1)\beta-(U+V)]}{(U-\beta)(V-\beta)}.
$$

The denominator is positive. The positive-beta return-before-send region is

$$
\boxed{\beta>\beta_c=\frac{U+V}{UV+1},\qquad
1-\beta_c=\frac{(U-1)(V-1)}{UV+1}.}
$$

For a proper delay, R<1 precisely when

$$
d<D(\beta)=\frac1B-A
=\frac{\beta[(UV+1)\beta-(U+V)]}{V\sqrt{1-\beta^2}(U-\beta)}.
$$

No d>=0 creates a loop outside the immediate region. Within that region there is a unique threshold for every finite d. To see this, for any speed q>1,

$$
\partial_\beta\left(\frac{q\sqrt{1-\beta^2}}{q-\beta}\right)
=\frac{q(1-q\beta)}{\sqrt{1-\beta^2}(q-\beta)^2}.
$$

Also beta_c>1/U and beta_c>1/V. Thus A and B strictly decrease beyond beta_c, so B(A+d) decreases from 1+d*B(beta_c) to zero as beta approaches 1. For d>0 the unique root beta_d lies strictly between beta_c and 1; for d=0 it is beta_c. A finite delay does not eliminate this region over the whole range beta<1.

This proof uses the decreasing branch. The return ratio is not globally decreasing in beta: R'(0)=1/U+(1+d)/V>0.

The reproduction code solves the root in e=1-beta rather than subtracting a nearly unit beta from 1:

$$
R(e)=\frac{UV e(2-e)}{(U-1+e)(V-1+e)}
+\frac{dV\sqrt{e(2-e)}}{V-1+e}.
$$

Its bracket is 0<=e<=(U-1)(V-1)/(UV+1), on which R(e) increases. Squaring the equation is unnecessary and could introduce extraneous roots. Decimal arithmetic preserves a stress-case threshold gap about 5e-31 that binary float cannot represent as 1-beta.

## Reply direction matters when a delay is present

The immediate result is symmetric in U and V. For a delayed reply,

$$
R(U,V)-R(V,U)=
\frac{d\beta\sqrt{1-\beta^2}(U-V)}{(U-\beta)(V-\beta)}.
$$

For a fixed unordered pair of speeds, assigning the faster speed to the reply lowers the arrival ratio when beta,d>0. This is a consequence of where the proper waiting interval is placed.

An exact example demonstrates a change of outcome:

| beta | d | U outbound | V return | t2/T | Reply before sending? |
|---|---|---|---|---|---|
| 4/5 | 1/5 | 3 | 2 | 56/55 | No |
| 4/5 | 1/5 | 2 | 3 | 54/55 | Yes |

Both immediate ratios are 9/11. The example was identified during the independent derivation, disclosed in a registration amendment and then confirmed with exact rational arithmetic.

## Light-speed boundaries and finite motion budgets

If the outbound leg is light-speed, U=1 and V>=1, the immediate ratio is V*(1+beta)/(V-beta)>=1. If the return leg is light-speed, V=1 and U>=1, it is U*(1+beta)/(U-beta)>=1. The nonnegative delay adds a nonnegative term. Therefore neither configuration yields a reply before sending in this construction. This is a statement about these two collinear inertial legs, not a universal result for every signaling geometry.

For an observer-speed cap beta_max<1 inside the immediate region, put A_cap=A(beta_max), B_cap=B(beta_max). The interval [0,beta_max] contains no return-before-send event exactly when

$$
\boxed{d\ge d_{\rm cap}=B_{\rm cap}^{-1}-A_{\rm cap}.}
$$

The decreasing-branch proof above establishes the interval statement. If beta_max<=beta_c, no delay is needed. Equality at the cap gives t2=T, not a strict loop. A fixed absolute tau does not uniformly enforce a fixed d when T is allowed to grow.

For a massive observer of fixed rest mass M, a kinetic-energy budget in Alice's frame supplies a separate standard-SR motion cap:

$$
\gamma_{\max}=1+\frac{E_{\rm kin,max}}{Mc^2},\qquad
\beta_{\max}=\sqrt{1-\gamma_{\max}^{-2}}.
$$

The modeled loop requires gamma>gamma_d=1/sqrt(1-beta_d^2). The JSON reports gamma_d and gamma_d-1 in units of observer rest energy. For U=V=2,d=100, gamma_d is about 200.0175, so this particular observer-motion threshold is about 199.0175*M*c^2. This excludes signal generation, propulsion, acceleration history, losses and apparatus costs; it is not an energy budget for physically producing an FTL signal.

## Near-light-speed and long-delay limits

For fixed U,V>1 and e=1-beta approaching zero,

$$
R=d\frac{\sqrt2 V}{V-1}\sqrt e
+\frac{2UV}{(U-1)(V-1)}e+O(d e^{3/2}+e^2).
$$

For fixed speeds and d tending to infinity, the threshold obeys

$$
\boxed{1-\beta_d\sim\frac{(V-1)^2}{2V^2d^2},\qquad
\gamma_d\sim\frac{dV}{V-1}.}
$$

The reply speed controls this leading long-delay term. The limit is not uniform as either speed approaches 1: the exact light-speed boundaries above have no loop region. The observer-motion cost diverges as the threshold approaches beta=1 for fixed positive M.

These calculations do not define a new observed population or speed bound. A population diagnostic for variable U,V,d would need its joint distribution specified. When both U and V tend to infinity at fixed d and a speed-independent beta population, the preceding delayed ceiling is recovered; finite unequal speeds cannot silently be substituted into the original one-parameter bound.

## A preferred time changes the operational model

Bob's return leg is future-directed in Bob's time, but its Alice-frame duration is

$$
\Delta t=\gamma(1-\beta V)\Delta t',\qquad\Delta t'>0.
$$

Every loop in this construction has beta>beta_c>1/V, making this return duration negative in Alice time. If Alice time were designated a global preferred time and every available signal and ordinary waiting operation had to increase it, this assumed return leg would be forbidden. A finite chain of strictly increasing preferred times cannot close.

This is an alternative signaling rule. It does not derive a preferred frame, Lorentz-invariant field dynamics, stability, positive energy, unitarity or a measurable front velocity. In curved spacetime or a theory with several effective propagation cones, a common global time compatible with every operational channel must actually exist; assigning separate local times does not prove that.

The primary connection is Liberati, Sonego and Visser, *Faster-than-c signals, special relativity, and causality*, [arXiv:gr-qc/0107091v2](https://arxiv.org/abs/gr-qc/0107091v2), version dated 14 February 2002, [DOI 10.1006/aphy.2002.6233](https://doi.org/10.1006/aphy.2002.6233). [Section 3.1](https://arxiv.org/html/gr-qc/0107091v2#S3.SS1) discusses the preferred-frame restriction; [section 3.2.2](https://arxiv.org/html/gr-qc/0107091v2#S3.SS2.SSS2) discusses a global temporal function and effective-metric stable causality. The exact-version abstract and HTML were live retrieved on 2026-10-02 UTC; [SOURCE_VERIFICATION.json](../research/asymmetric-delay/SOURCE_VERIFICATION.json) records URLs, timestamps and hashes. This was a targeted primary-source check, not an exhaustive or current-discovery literature survey. No article text or private archive file is redistributed here.

## Reproduction and limits

From the repository root:

```bash
python3 research/asymmetric-delay/reproduce.py --out /tmp/ftl-asymmetric-delay-results.json
```

Python >=3.10 and only the standard library are required; Python 3.12.14 was used. The [registration](../research/asymmetric-delay/REGISTRATION.md) and its disclosed amendment fix checks and tolerances. [RESULTS.json](../research/asymmetric-delay/RESULTS.json) records parameters, hash, residuals and failure status. The script exits nonzero on a failed check.

The run passed 5,000 direct event comparisons, 288 exact rational events, 128 high-precision thresholds, 256 side checks, 25 diagonal comparisons, 25 caps with 6,425 grid points, 150 light-speed boundary cases, nine long-delay speed pairs, and 128 preferred-time signs and observer-energy identities. Numerical checks complement the threshold and cap proofs. See the separate [internal audit](../research/asymmetric-delay/AUDIT.md) for independent routes and limitations.

Energy and stability of a hypothetical propagation mechanism remain unresolved. So do a derived FTL-to-PTA/ringdown observable, authenticated posterior inputs and calibrated diagnostic interpretation. The equations characterize an assumed controllable spacelike signaling rule; they do not show nature implements it.
