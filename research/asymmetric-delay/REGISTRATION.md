# Unequal signal speeds and proper reply delay: replication protocol

Recorded 2026-10-01 Pacific (2026-10-02 UTC), before the accompanying decisive numerical run. This is subsequent conditional mathematical work; the historical snapshot remains unchanged.

Prior knowledge: the equal-speed delay model was reproduced in the preceding round. Algebra for the unequal-speed candidate was considered before this registration, and a separate internal derivation was requested. This is a disclosed replication and sensitivity protocol, not a blind preregistration or novelty claim. All inputs are generated mathematical parameters; none are observations or private raw files.

Model: collinear inertial observers in Minkowski spacetime, coincident origins, Alice sending at T>0. Outbound speed U*c is controllable in Alice's frame; return speed V*c is controllable in Bob's frame. Initially U,V>1, 0<=beta<1, and Bob's proper wait is tau>=0 with d=tau/T. Boundary tests allow U=1 or V=1. No energy, production, stability, preferred-frame dynamics or empirical model is supplied.

Candidate: s=sqrt(1-beta^2), A=U*s/(U-beta), B=V*s/(V-beta), t2/T=B*(A+d). Immediate threshold beta_c=(U+V)/(U*V+1); its stable gap is (U-1)*(V-1)/(U*V+1). Delayed thresholds will be solved in the gap e=1-beta, on [0,1-beta_c], without squaring the root equation. No global monotonicity in beta is presumed.

Fixed checks and failure criteria:

1. Compare 5,000 direct Lorentz-event calculations with B*(A+d), seed 20261002; U,V=1+10^z, z uniform on [-3,2], beta uniform on [0,0.999], d=0 every fifth case and otherwise 10^z with z uniform on [-3,2]. Require absolute/relative agreement 1e-10. Also compare exact Fraction events for six rational beta/s pairs, four U values, four V values and three delays (288 cases).
2. Verify the exact example U=2,V=3,beta=12/13,d=1/2 gives 85/126, and swapping the speeds gives 95/126. Immediate arrival must be symmetric under U/V exchange; the delayed term need not be. Check recovery of the preceding equal-speed formula and threshold.
3. Use 90-digit Decimal arithmetic for U,V in {1.001,1.1,2,10,100}, d in {0,0.1,1,10,100}: 125 thresholds. Require positive gaps within the immediate gap, root residual below 1e-60 and the correct strict signs at 0.999 and 1.001 times the delayed gap, whenever inside [0,1]. For the 25 diagonal cases compare with the earlier exact equal-speed gap, relative tolerance 1e-60. Add three narrow-gap stress cases, including U=1.000000000001,V=1.000000001,d=1000000. Decimal tests must retain positive gaps rather than subtracting beta from 1 in binary float.
4. For each of the 25 speed pairs, put the cap halfway into the immediate region in gap coordinates. Verify d_cap=1/B_cap-A_cap gives t2/T=1 to 1e-60; 0.999*d_cap permits a loop at the cap, whereas 1.001*d_cap excludes it at 257 uniformly spaced beta values through the cap, tolerance 1e-12. The interval conclusion also requires a mathematical monotonicity proof on the immediate region; grid checks alone are insufficient.
5. Test light-speed boundary U=1 or V=1 with the other speed in {1,1.1,2,10,100}, beta in {0,0.1,0.5,0.9,0.999}, d in {0,0.1,10}. Require t2/T>=1 within 1e-12. Do not extrapolate this statement to different geometries or general theories.
6. For fixed U,V in {1.1,2,10}, test the large-d gap coefficient d^2*(1-beta_d) -> (V-1)^2/(2*V^2), using d=1000 and 1000000. Require the latter relative error below 1e-6 and closer to the limit than the former. This is a fixed-U,V limit; no interchange with U approaching 1 is assumed.
7. A separate internal audit must check events, threshold uniqueness, speed-cap conditions, asymptotic assumptions and physical claim limits. In the loop region, test the sign 1-beta*V<0 for the return leg's Alice-frame duration; explain the consequences if Alice time were chosen as a globally increasing preferred time. That optional rule is a changed operational model, not derived physical dynamics.

Record Python version, script hash, parameters, counts, residuals, tolerances and all failures in JSON. A failed or skipped decisive check is not success; the script must return nonzero on any failure. Python >=3.10, standard library only. Reproduce with `python3 research/asymmetric-delay/reproduce.py --out /tmp/ftl-asymmetric-delay-results.json` from the repository root. Preserve this registration; append dated amendments instead of rewriting it after the run.

## Amendment 1 — 2026-10-02 UTC, after the initial run

The initial registered checks passed. A separate internal derivation then supplied a sharper exact speed-exchange example, already known before this amendment: beta=4/5,d=1/5 gives R(3,2)=56/55>1 and R(2,3)=54/55<1. Add this exact Fraction check to the next run, preserving the original protocol and check counts. This is a disclosed counterexample confirmation.

Also report the standard massive-observer Lorentz threshold gamma_d=1/sqrt(e_d*(2-e_d)) and kinetic energy in units of M*c^2, gamma_d-1. Cross-check gamma_d^2*(1-beta_d^2)=1 at the existing 128 Decimal thresholds, tolerance 1e-60. This adds a conditional observer-motion cost; it is not an energy estimate for a signal-producing mechanism, and does not imply an attainable FTL channel. The definition was considered before registering this added check.
