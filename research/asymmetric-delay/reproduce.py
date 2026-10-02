#!/usr/bin/env python3
"""Unequal hypothetical signal speeds with Bob's proper reply delay.

Python >=3.10, standard library only; mathematical parameters, no observations.
See REGISTRATION.md and ../../docs/ASYMMETRIC_SIGNAL_SPEEDS.md.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import random


SPEEDS = ("1.001", "1.1", "2", "10", "100")
DELAYS = ("0", "0.1", "1", "10", "100")
STRESS = (("1.000000000001", "1.000000001", "1000000"),
          ("1.000000001", "100", "1000000"),
          ("100", "1.000000001", "1000000"))


def factor(speed: float, beta: float) -> float:
    return speed * math.sqrt((1 - beta) * (1 + beta)) / (speed - beta)


def return_ratio(u: float, v: float, beta: float, delay: float) -> float:
    return factor(v, beta) * (factor(u, beta) + delay)


def event_ratio(u: float, v: float, beta: float, delay: float) -> float:
    # T=c=1. Use four event transforms/intersections rather than the factors.
    gamma = 1 / math.sqrt(1 - beta * beta)
    reception_t = u / (u - beta)
    reception_x = beta * reception_t
    reception_bob_t = gamma * (reception_t - beta * reception_x)
    departure_bob_t = reception_bob_t + delay
    arrival_bob_t = v * departure_bob_t / (v - beta)
    arrival_bob_x = -beta * arrival_bob_t
    return gamma * (arrival_bob_t + beta * arrival_bob_x)


def ratio_at_gap(u: Decimal, v: Decimal, delay: Decimal, gap: Decimal) -> Decimal:
    s = (gap * (2 - gap)).sqrt()
    a = u * s / (u - 1 + gap)
    b = v * s / (v - 1 + gap)
    return b * (a + delay)


def threshold_gap(u: Decimal, v: Decimal, delay: Decimal) -> Decimal:
    """Root on the proved monotone loop branch, with no subtraction of beta."""
    upper = (u - 1) * (v - 1) / (u * v + 1)
    if not (u > 1 and v > 1 and delay >= 0):
        raise ValueError("Threshold requires U,V>1 and d>=0")
    if not delay:
        return upper
    lower = Decimal(0)
    for _ in range(320):
        middle = (lower + upper) / 2
        if middle == lower or middle == upper:
            break
        if ratio_at_gap(u, v, delay, middle) < 1:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def equal_speed_gap(u: Decimal, delay: Decimal) -> Decimal:
    # Closed-form diagonal reference from the preceding round.
    h = 2 / ((delay * delay + 4).sqrt() + delay)
    radical = (u * u * (1 - h * h) + h * h).sqrt()
    return h * h * (u - 1) ** 2 / (u * u - h * h * (u - 1) + u * radical)


def run_checks() -> dict:
    failures = []
    rng = random.Random(20261002)
    event_errors, swap_errors, diagonal_errors = [], [], []
    for i in range(5000):
        u, v = (1 + 10 ** rng.uniform(-3, 2) for _ in range(2))
        beta = rng.uniform(0, .999)
        delay = 0.0 if i % 5 == 0 else 10 ** rng.uniform(-3, 2)
        actual, expected = event_ratio(u, v, beta, delay), return_ratio(u, v, beta, delay)
        event_errors.append(abs(actual - expected) / expected)
        if not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10):
            failures.append(f"Lorentz event mismatch {i}")
        difference = expected - return_ratio(v, u, beta, delay)
        analytic = delay * math.sqrt((1 - beta) * (1 + beta)) * beta * (u - v) / ((u - beta) * (v - beta))
        swap_errors.append(abs(difference - analytic) / max(1, abs(expected)))
        if not math.isclose(difference, analytic, rel_tol=1e-10, abs_tol=1e-10):
            failures.append(f"Speed-exchange identity {i}")
        a = factor(u, beta)
        diagonal_errors.append(abs(return_ratio(u, u, beta, delay) - (a * a + delay * a)))
        if not math.isclose(return_ratio(u, u, beta, delay), a * a + delay * a, rel_tol=1e-10, abs_tol=1e-10):
            failures.append(f"Equal-speed recovery {i}")

    exact_cases = 0
    for beta, s in ((Fraction(3, 5), Fraction(4, 5)), (Fraction(4, 5), Fraction(3, 5)),
                    (Fraction(5, 13), Fraction(12, 13)), (Fraction(12, 13), Fraction(5, 13)),
                    (Fraction(8, 17), Fraction(15, 17)), (Fraction(15, 17), Fraction(8, 17))):
        for u in map(Fraction, ("3/2", "2", "3", "10")):
            for v in map(Fraction, ("3/2", "2", "3", "10")):
                for delay in map(Fraction, ("0", "1/2", "2")):
                    reception_t = u / (u - beta)
                    reception_bob_t = (reception_t - beta * beta * reception_t) / s
                    arrival_bob_t = v * (reception_bob_t + delay) / (v - beta)
                    actual = (arrival_bob_t - beta * beta * arrival_bob_t) / s
                    a, b = u * s / (u - beta), v * s / (v - beta)
                    exact_cases += 1
                    if actual != b * (a + delay):
                        failures.append("Fraction event mismatch")
    a, b, delay = Fraction(5, 7), Fraction(5, 9), Fraction(1, 2)
    if b * (a + delay) != Fraction(85, 126) or a * (b + delay) != Fraction(95, 126):
        failures.append("Exact speed-exchange example")
    beta, s, d = Fraction(4, 5), Fraction(3, 5), Fraction(1, 5)
    a, b = 3 * s / (3 - beta), 2 * s / (2 - beta)
    if b * (a + d) != Fraction(56, 55) or a * (b + d) != Fraction(54, 55):
        failures.append("Exact loop/no-loop speed-exchange counterexample")

    rows, caps, asymptotics = [], [], []
    side_cases, diagonal_cases, cap_grid_cases, preferred_sign_cases, energy_cases = 0, 0, 0, 0, 0
    max_residual, max_equal_gap_error, max_cap_residual = Decimal(0), Decimal(0), Decimal(0)
    with localcontext() as context:
        context.prec = 90
        cases = [(u, v, d) for u in SPEEDS for v in SPEEDS for d in DELAYS] + list(STRESS)
        for u_text, v_text, d_text in cases:
            u, v, d = map(Decimal, (u_text, v_text, d_text))
            gap = threshold_gap(u, v, d)
            immediate_gap = (u - 1) * (v - 1) / (u * v + 1)
            residual = abs(ratio_at_gap(u, v, d, gap) - 1)
            max_residual = max(max_residual, residual)
            if not 0 < gap <= immediate_gap or residual > Decimal("1e-60"):
                failures.append(f"Threshold root {u_text},{v_text},{d_text}")
            for multiplier, precedes in ((Decimal(".999"), True), (Decimal("1.001"), False)):
                side_cases += 1
                if (ratio_at_gap(u, v, d, gap * multiplier) < 1) != precedes:
                    failures.append(f"Threshold side {u_text},{v_text},{d_text}")
            if u == v:
                diagonal_cases += 1
                relative = abs(gap - equal_speed_gap(u, d)) / gap
                max_equal_gap_error = max(max_equal_gap_error, relative)
                if relative > Decimal("1e-60"):
                    failures.append(f"Equal-speed threshold {u_text},{d_text}")
            preferred_sign_cases += 1
            beta_inside = 1 - gap / 2
            if 1 - beta_inside * v >= 0:
                failures.append(f"Return leg preferred-time sign {u_text},{v_text},{d_text}")
            gamma = 1 / (gap * (2 - gap)).sqrt()
            energy_cases += 1
            if abs(gamma * gamma * gap * (2 - gap) - 1) > Decimal("1e-60"):
                failures.append(f"Observer Lorentz factor {u_text},{v_text},{d_text}")
            rows.append({"U": u_text, "V": v_text, "d": d_text,
                         "gap_1_minus_beta_threshold": str(gap), "root_residual": str(residual),
                         "gamma_threshold": str(gamma), "observer_kinetic_energy_over_Mc2": str(gamma - 1)})

        for u_text in SPEEDS:
            for v_text in SPEEDS:
                u, v = Decimal(u_text), Decimal(v_text)
                gap = (u - 1) * (v - 1) / (u * v + 1) / 2
                s = (gap * (2 - gap)).sqrt()
                a, b = u * s / (u - 1 + gap), v * s / (v - 1 + gap)
                d_cap = 1 / b - a
                residual = abs(ratio_at_gap(u, v, d_cap, gap) - 1)
                max_cap_residual = max(max_cap_residual, residual)
                if d_cap <= 0 or residual > Decimal("1e-60"):
                    failures.append(f"Cap equality {u_text},{v_text}")
                if ratio_at_gap(u, v, Decimal(".999") * d_cap, gap) >= 1:
                    failures.append(f"Cap lower-delay sign {u_text},{v_text}")
                cap = float(1 - gap)
                for i in range(257):
                    cap_grid_cases += 1
                    beta = cap * i / 256
                    if return_ratio(float(u), float(v), beta, float(Decimal("1.001") * d_cap)) < 1 - 1e-12:
                        failures.append(f"Cap grid {u_text},{v_text},{i}")
                caps.append({"U": u_text, "V": v_text, "cap_gap_1_minus_beta": str(gap),
                             "d_cap": str(d_cap), "equality_residual": str(residual)})

        for u_text in ("1.1", "2", "10"):
            for v_text in ("1.1", "2", "10"):
                u, v = Decimal(u_text), Decimal(v_text)
                limit = (v - 1) ** 2 / (2 * v * v)
                values = []
                for d_text in ("1000", "1000000"):
                    d = Decimal(d_text)
                    value = d * d * threshold_gap(u, v, d)
                    values.append({"d": d_text, "d2_times_gap": str(value),
                                   "relative_error": str(abs(value - limit) / limit)})
                if Decimal(values[1]["relative_error"]) > Decimal("1e-6") or not Decimal(values[1]["relative_error"]) < Decimal(values[0]["relative_error"]):
                    failures.append(f"Large-delay limit {u_text},{v_text}")
                asymptotics.append({"U": u_text, "V": v_text, "limit": str(limit), "values": values})

    boundary_cases = 0
    min_boundary_ratio = math.inf
    for other in (1.0, 1.1, 2.0, 10.0, 100.0):
        for u, v in ((1.0, other), (other, 1.0)):
            for beta in (0.0, .1, .5, .9, .999):
                for d in (0.0, .1, 10.0):
                    boundary_cases += 1
                    ratio = return_ratio(u, v, beta, d)
                    min_boundary_ratio = min(min_boundary_ratio, ratio)
                    if ratio < 1 - 1e-12:
                        failures.append(f"Light-speed boundary {u},{v},{beta},{d}")

    counts = {"floating_event_cases": len(event_errors), "exact_event_cases": exact_cases,
              "threshold_cases": len(rows), "threshold_side_cases": side_cases,
              "equal_speed_threshold_cases": diagonal_cases, "cap_cases": len(caps),
              "cap_grid_cases": cap_grid_cases, "light_speed_boundary_cases": boundary_cases,
              "large_delay_speed_pairs": len(asymptotics), "preferred_time_sign_cases": preferred_sign_cases,
              "observer_energy_cases": energy_cases}
    if counts != {"floating_event_cases": 5000, "exact_event_cases": 288, "threshold_cases": 128,
                  "threshold_side_cases": 256, "equal_speed_threshold_cases": 25, "cap_cases": 25,
                  "cap_grid_cases": 6425, "light_speed_boundary_cases": 150,
                  "large_delay_speed_pairs": 9, "preferred_time_sign_cases": 128, "observer_energy_cases": 128}:
        failures.append("Registered check counts not met")
    return {"run_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "claim_status": "conditional kinematics; no novelty, physical FTL or observational claim",
            "inputs": "deterministic mathematical parameters; no raw research data",
            "parameters": {"seed": 20261002, "threshold_U_and_V": list(SPEEDS), "threshold_delays": list(DELAYS),
                           "stress_U_V_d": list(STRESS), "decimal_precision": 90, "bisection_iterations_max": 320,
                           "cap_grid_points": 257},
            "tolerances": {"event_relative_and_absolute": 1e-10, "threshold_and_cap_residual": "1e-60",
                           "diagonal_gap_relative": "1e-60", "boundary_and_cap_grid_absolute": 1e-12,
                           "large_delay_final_relative": "1e-6"},
            "passed": not failures, "failures": failures, "checks": counts,
            "residuals": {"max_event_relative_error": max(event_errors), "max_speed_exchange_scaled_error": max(swap_errors),
                          "max_equal_speed_ratio_absolute_error": max(diagonal_errors), "max_threshold_residual": str(max_residual),
                          "max_equal_speed_gap_relative_error": str(max_equal_gap_error), "max_cap_residual": str(max_cap_residual),
                          "min_light_speed_boundary_ratio": min_boundary_ratio},
            "thresholds": rows, "speed_caps": caps, "large_delay_gap": asymptotics}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run_checks()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: result[k] for k in ("passed", "failures", "checks", "residuals")}))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
