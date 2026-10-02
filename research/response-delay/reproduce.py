#!/usr/bin/env python3
"""Finite proper reply delay in a hypothetical SR signaling rule.

Python >=3.10; standard library only. Generated parameters are not observations.
See REGISTRATION.md and ../../docs/FINITE_RESPONSE_DELAY.md for assumptions.
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


def factor(w: float, beta: float) -> float:
    return w * math.sqrt((1 - beta) * (1 + beta)) / (w - beta)


def return_ratio(w: float, beta: float, delay: float) -> float:
    a = factor(w, beta)
    return a * a + delay * a


def event_ratio(w: float, beta: float, delay: float) -> float:
    # Intersect the outgoing ray with Bob, transform that reception event,
    # wait in Bob's frame, intersect the return ray with Alice, and transform.
    gamma = 1 / math.sqrt(1 - beta * beta)
    reception_t = w / (w - beta)
    reception_x = beta * reception_t
    reception_bob_t = gamma * (reception_t - beta * reception_x)
    departure_bob_t = reception_bob_t + delay
    arrival_bob_t = w * departure_bob_t / (w - beta)
    arrival_bob_x = -beta * arrival_bob_t
    return gamma * (arrival_bob_t + beta * arrival_bob_x)


def threshold_gap(w, delay):
    """Stable exact expression for 1-beta_threshold; accepts float or Decimal."""
    if isinstance(w, Decimal):
        one = Decimal(1)
        h = 2 / ((delay * delay + 4).sqrt() + delay)
        radical = (w * w * (one - h * h) + h * h).sqrt()
    else:
        one = 1.0
        h = 2 / (math.hypot(delay, 2) + delay)
        radical = math.sqrt(w * w * (one - h * h) + h * h)
    return h * h * (w - one) ** 2 / (w * w - h * h * (w - one) + w * radical)


def simpson(function, end: float, n: int) -> float:
    values = [(1 if i in (0, n) else 4 if i % 2 else 2) * function(end * i / n)
              for i in range(n + 1)]
    return end * math.fsum(values) / (3 * n)


def ceiling(delay: float, n: int = 8192) -> float:
    h = 2 / (math.hypot(delay, 2) + delay)
    end = math.asin(h)
    return simpson(lambda u: (math.cos(u) ** 2 - delay * math.sin(u)) ** 2 * math.sin(u), end, n)


def severity(w: float, delay: float, n: int = 8192) -> float:
    gap = threshold_gap(w, delay)
    # acos(1-gap) would lose precision near beta=1.
    end = 2 * math.asin(math.sqrt(gap / 2))

    def integrand(u):
        # w-cos(u) is evaluated without subtracting two nearly equal numbers.
        a = w * math.sin(u) / ((w - 1) + 2 * math.sin(u / 2) ** 2)
        positive = max(1 - a * a - delay * a, 0.0)
        return positive * positive * math.sin(u)

    return simpson(integrand, end, n)


def run_checks() -> dict:
    failures = []
    rng = random.Random(20261001)
    relative_errors = []
    for i in range(5000):
        w = 1 + 10 ** rng.uniform(-3, 2)
        beta = rng.uniform(0, .999)
        delay = 0.0 if i % 5 == 0 else 10 ** rng.uniform(-3, 2)
        actual, expected = event_ratio(w, beta, delay), return_ratio(w, beta, delay)
        relative_errors.append(abs(actual - expected) / expected)
        if not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10):
            failures.append(f"Lorentz event mismatch at case {i}")

    exact_cases = 0
    for beta, s in [(Fraction(3, 5), Fraction(4, 5)), (Fraction(4, 5), Fraction(3, 5)),
                    (Fraction(5, 13), Fraction(12, 13)), (Fraction(12, 13), Fraction(5, 13)),
                    (Fraction(8, 17), Fraction(15, 17)), (Fraction(15, 17), Fraction(8, 17))]:
        for w in map(Fraction, ("3/2", "2", "3", "10")):
            for delay in map(Fraction, ("0", "1/2", "2")):
                reception_t = w / (w - beta)
                reception_x = beta * reception_t
                reception_bob_t = (reception_t - beta * reception_x) / s
                arrival_bob_t = w * (reception_bob_t + delay) / (w - beta)
                event = (arrival_bob_t - beta * beta * arrival_bob_t) / s
                a = w * s / (w - beta)
                exact_cases += 1
                if event != a * a + delay * a:
                    failures.append("Exact Fraction event mismatch")

    a = Fraction(5, 7)
    if a * a + Fraction(1, 2) * a != Fraction(85, 98) or 1 / a - a != Fraction(24, 35):
        failures.append("Exact documented example failed")

    threshold_cases, side_cases, max_residual = 0, 0, Decimal(0)
    with localcontext() as context:
        context.prec = 70
        for w_text in ("1.001", "1.1", "2", "10", "100"):
            for d_text in ("0", "0.1", "0.5", "1", "10", "100"):
                w, delay = Decimal(w_text), Decimal(d_text)
                gap = threshold_gap(w, delay)

                def ratio_at_gap(e):
                    a = w * (e * (2 - e)).sqrt() / (w - 1 + e)
                    return a * a + delay * a

                residual = abs(ratio_at_gap(gap) - 1)
                max_residual = max(max_residual, residual)
                threshold_cases += 1
                if residual > Decimal("1e-45"):
                    failures.append(f"Threshold residual W={w}, d={delay}")
                for multiplier, before_send in ((Decimal(".999"), True), (Decimal("1.001"), False)):
                    side_cases += 1
                    if (ratio_at_gap(gap * multiplier) < 1) != before_send:
                        failures.append(f"Threshold side W={w}, d={delay}")

    ceiling_rows = []
    max_resolution_change = 0.0
    for delay in (0.0, .1, .5, 1.0, 2.0, 10.0, 100.0, 1000.0):
        low, high = ceiling(delay, 4096), ceiling(delay, 8192)
        change = abs(high - low) / high
        max_resolution_change = max(max_resolution_change, change)
        if change > 1e-9 or high <= 0:
            failures.append(f"Ceiling quadrature d={delay}")
        ceiling_rows.append({"delay_tau_over_T": delay, "J_infinity_uniform_beta": high,
                             "delay2_times_J": delay * delay * high})
    if abs(ceiling_rows[0]["J_infinity_uniform_beta"] - .2) > 1e-12:
        failures.append("Immediate ceiling failed")
    if abs(ceiling_rows[-1]["delay2_times_J"] - 1 / 12) > 1e-6:
        failures.append("Large-delay asymptotic failed")
    if not all(a["J_infinity_uniform_beta"] > b["J_infinity_uniform_beta"]
               for a, b in zip(ceiling_rows, ceiling_rows[1:])):
        failures.append("Ceiling must strictly decrease with delay")

    finite_rows = []
    for delay in (0.0, .1, .5, 1.0, 2.0, 10.0):
        limit = ceiling(delay)
        values = [severity(w, delay) for w in (1.1, 2, 10, 100, 1e6)]
        if not all(0 < value <= limit * (1 + 1e-12) for value in values):
            failures.append(f"Finite severity above ceiling or nonpositive d={delay}")
        if not all(a < b for a, b in zip(values, values[1:])):
            failures.append(f"Finite severity not increasing with W d={delay}")
        if abs(values[-1] - limit) / limit > 1e-4:
            failures.append(f"Finite severity not approaching ceiling d={delay}")
        finite_rows.append({"delay_tau_over_T": delay, "W": [1.1, 2, 10, 100, 1e6], "J": values})

    return {
        "run_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "claim_status": "conditional mathematical sensitivity; no observations or physical FTL claim",
        "inputs": "deterministic generated parameters; fixed uniform-beta population and fixed tau/T",
        "passed": not failures, "failures": failures,
        "checks": {"floating_event_cases": 5000, "max_event_relative_error": max(relative_errors),
                   "exact_event_cases": exact_cases, "decimal_threshold_cases": threshold_cases,
                   "threshold_side_cases": side_cases, "max_threshold_residual": str(max_residual),
                   "ceiling_resolutions": [4096, 8192], "max_ceiling_resolution_change": max_resolution_change},
        "uniform_ceiling": ceiling_rows, "finite_W_uniform_severity": finite_rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run_checks()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"passed": result["passed"], "failures": result["failures"], "checks": result["checks"]}))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
