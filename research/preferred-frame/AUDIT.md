# Internal audit of the fixed-preferred-time scalar exercise

2026-10-01 Pacific (2026-10-02 UTC). Separate mathematical and source reviews checked a newly assumed free-field toy model. This is internal AI-assisted checking, not external peer review, blind discovery preregistration, empirical replication or a novelty claim.

The final [root script](reproduce.py) and [result](RESULTS.json) passed the [registered protocol](REGISTRATION.md) and its disclosed amendments. Amendment1 corrected a finite-difference scale before any decisive run. Amendment2 records that initial literal pathology/discriminant outputs and frequency-only energy checking did not constitute sufficient control coverage; the final script evaluates shared signed dispersion, Hamiltonian and principal formulas, compares healthy controls, and gates dispersion/average-energy residuals. It also corrects the initial (-,+,+,+) covector-norm wording. Only the current hash-bound result supports these strengthened checks. No preferred-model physical instability was discovered or hidden.

Root validation: 2,000 deterministic modes, 6,000 dispersion derivative components, 2,000 boosted on-shell modes, 21 exact principal matrices, 12 compact-wave energy cells at 256/512 points, 2,000 preferred relay constructions, 20 laboratory timing cells, four small-velocity and five light-speed boundaries. All pathological controls were correctly classified. Maximum group-derivative absolute error was 2.84e-10; boosted invariant scaled error 3.35e-12; compact-wave energy relative error 3.13e-11; worldline/timing relative disagreements below 4e-15. Parameters, tolerances, Python version and SHA are saved. These finite checks supplement the analytic arguments; they are not a theorem about an unknown interacting theory.

The separately formulated [independent script](independent_check.py) uses exact rational Pythagorean boosts and direct event intersections. Its [result](INDEPENDENT_RESULTS.json) passes 9,206 checks of transformed principal symbols, on-shell frequencies, discriminant squares, preferred translation-charge square completion, finite-propagation flux bounds, longitudinal intersections, transverse squared intersections and orientation signs. It imports no root implementation. The relocated code adds a run timestamp; its mathematical algorithms are unchanged. Python 3.12.14 was used; Python >=3.10, standard library only. Run with `python3 research/preferred-frame/independent_check.py --out /tmp/ftl-preferred-frame-independent.json` from the public repository root.

The strongest claim distinctions were checked explicitly:

- Free preferred H>=0 and real frequencies hold for positive kinetic/gradient signs and m^2>=0. Arbitrary prescribed forcing or completed interactions do not inherit a supplied stability or source-energy proof.
- Constant preferred time orders every permitted faster-cone signal and ordinary relay. This sufficient no-loop condition changes the original reciprocal sender-frame assumption. The loop diagnostic consequently vanishes in this prescription; no universal speed limit or physical observation follows.
- Front=a holds for the complete local continuum PDE. A low-energy action with unspecified ultraviolet completion cannot settle the physical front.
- t'=constant ceases to be a Cauchy surface at |v|*a>=1. Arbitrary-data transverse growth there is a wrong-slicing obstruction. Positive transformed physical frequencies and translation charges on admissible data are distinct from signed formal canonical flux on non-Cauchy surfaces.
- Ordinary rods/clocks, controllable massless pulses, source/detector coupling and latency calibration are explicit laboratory assumptions. The calculated roundtrip anisotropy is a conditional prediction, not a measurement or usable apparatus.

Primary-source review verified exact-version arXiv texts and metadata for Liberati–Sonego–Visser, Bruneton and Babichev–Mukhanov–Vikman. See [SOURCE_VERIFICATION.json](SOURCE_VERIFICATION.json). The papers' background, finite-clump, conjecture and effective-theory assumptions are distinguished in the new note; no source establishes this toy model's physical viability. Cached article contents are excluded from published files.

The final skeptical check also uses a generic sphere/moving-target root solver rather than the closed-form lab durations. [GEOMETRY_RESULTS.json](GEOMETRY_RESULTS.json) records 2,400 preferred relays, 2,400 distinct sender-frame events, 1,600 arbitrary-angle laboratory geometries, 25 registered/light-speed cells and 21 exact tensor transformations. [FAILURE_GATE_RESULTS.json](FAILURE_GATE_RESULTS.json) records four in-memory defects rejected by the strengthened root checks: erased negative dispersion, erased negative energy, erased negative boosted coefficient and a perturbed frequency. Candidate files remain unchanged. These bounded checks test implemented failure gates, not all possible mutations or physical stability.

Reproduce the original reviewer algorithms (with bytecode writing disabled in the public wrappers):

```bash
python3 research/preferred-frame/geometry_check.py --candidate research/preferred-frame/reproduce.py --out /tmp/ftl-preferred-geometry.json
python3 research/preferred-frame/failure_gate_check.py --candidate research/preferred-frame/reproduce.py --out /tmp/ftl-preferred-gates.json
```

The private v6.3 kit/proofs, authenticated compatible posterior inputs, diagnostic calibration and physical observational bridge remain unavailable. No unrelated project, raw private archive, website deployment, archival DOI, strict empirical v1.3.0 result or continuously running agent is established by this round.
