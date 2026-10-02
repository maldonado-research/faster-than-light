# Internal audit of unequal-speed proper-delay kinematics

Recorded 2026-10-01 Pacific (2026-10-02 UTC). A separate internal derivation and skeptical workflow/mathematical review checked this methods extension. These are internal AI-assisted checks, not external peer review or independent physical replication. No observational data or private raw research files enter the calculation, and no mathematical novelty is claimed.

The [registered root run](REGISTRATION.md) and disclosed amendment passed the implemented checks in [RESULTS.json](RESULTS.json). Python 3.12.14 was used; code requires Python >=3.10 and only the standard library. Root results identify the current script hash, parameters, tolerances and failures. A second registered run added only the disclosed exact orientation counterexample and observer Lorentz-factor reporting.

## Different event and threshold routes

The root computes the outgoing intersection, transforms reception to Bob, waits in Bob time, intersects the reply with Alice in Bob coordinates, and transforms back. The independent route instead places Bob's delayed departure directly at t_r=t1+gamma*tau in Alice coordinates. It transforms the return displacement and intersects Alice there. The routes agree on R=B*(A+d).

The root solves the threshold using e=1-beta. The independent implementation uses q=sqrt[(1-beta)/(1+beta)], with A_S=2*S*q/[(S-1)+(S+1)*q^2]. Multiplying by positive denominators, without squaring, gives

$$
(U+1)(V+1)q^4-2dV(U+1)q^3-2(UV+1)q^2
-2dV(U-1)q+(U-1)(V-1)=0.
$$

It solves the unique loop-branch root between 0 and q_c=sqrt[(U-1)(V-1)/((U+1)(V+1))], then maps e=2*q^2/(1+q^2). At 110-digit precision, all 128 root gaps agree within 2.39e-76 relative, including the three stress cases. This is an independently formulated route to the same conditional model.

The derivative proof was checked separately: beta_c>1/U and beta_c>1/V, so both factors decrease throughout the immediate loop region. R increases with d and decreases with beta on that region, establishing the unique delayed threshold. Outside that region, R_0>=1 and d*B>=0. Thus the finite-cap conclusion follows algebraically; a grid alone would not prove it. R'(0)>0 rules out global beta monotonicity.

## Reproducible independent checks

From the repository root:

```bash
python3 research/asymmetric-delay/independent_check.py --out /tmp/ftl-asymmetric-independent.json --root-result research/asymmetric-delay/RESULTS.json
```

Read the fresh output and its exit status. The relocated audit script retains the independent mathematical algorithms; its default path and output provenance were made portable for this public record. [INDEPENDENT_RESULTS.json](INDEPENDENT_RESULTS.json) records its hash, root-result hash, root-script hash consistency and the comparisons.

Its executed checks passed:

- 20,000 Alice-coordinate event comparisons, seed 20261002, absolute/relative tolerance 2e-10; maximum relative disagreement 2.54e-13.
- 256 exact rational event and delay-allowance cases, plus the 56/55 versus 54/55 orientation counterexample.
- 150 Decimal thresholds and 300 side checks, with residuals below 1e-60 (observed maximum 1e-84 under 85-digit residual evaluation).
- 75,225 cap equality, side and grid checks, grid tolerance 2e-10.
- 120 light-speed-leg cases and 27 near-beta-one expansion cases.
- 12 long-delay comparisons over four unequal speed pairs; the final d=10000 discrepancies satisfy relative tolerance 1e-6.
- 128 comparisons to the stored root thresholds through the unsquared rapidity polynomial, relative tolerance 1e-60, using 110-digit arithmetic.

The independent test plan was developed after the initial algebra and in conversation with the root protocol. It was not a blind or externally preregistered audit. Exact roots, inequalities and limits are mathematical deductions within the assumptions; finite checks do not prove validity outside that model.

## Boundaries and source checking

Either light-speed leg removes the return-before-send region. The beta-to-1 expansion is nonuniform as U or V approaches 1. At fixed strictly superluminal U,V and large d, the gap coefficient is (V-1)^2/(2*V^2), and the observer Lorentz factor grows like d*V/(V-1). The reported kinetic energy uses a fixed positive observer rest mass; it supplies no signal-generation energy model.

In the loop region, the Bob-future return segment has Alice-coordinate duration proportional to 1-beta*V<0. Requiring every operational channel to increase a common global preferred time would reject it. This changes the available signaling rule. A coordinate label alone does not supply such dynamics or establish stable causality for several effective cones.

The preferred-time interpretation was compared with live retrieved Liberati–Sonego–Visser arXiv:gr-qc/0107091v2, sections 3.1 and 3.2.2. Retrieval metadata is in [SOURCE_VERIFICATION.json](SOURCE_VERIFICATION.json); article contents are not redistributed. This is an established primary-source connection, not an exhaustive search or evidence for a realizable FTL channel.

No failure occurred in these registered and independent runs. The earlier archived narrow-wedge failure remains recorded in its own checkpoint; these new calculations do not validate the missing private v6.3 toolchain. Physical fields, production, stability, unitarity, measurable predictions and authenticated observational inference remain unresolved.
