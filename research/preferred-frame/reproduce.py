#!/usr/bin/env python3
"""Conditional preferred-frame free scalar checks; Python >=3.10, stdlib only.

All inputs are deterministic mathematical parameters, not observations.
See REGISTRATION.md and ../../docs/PREFERRED_FRAME_SCALAR.md.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import random


def direction(rng):
    z, theta = rng.uniform(-1, 1), rng.uniform(0, 2 * math.pi)
    radial = math.sqrt(1 - z * z)
    return (radial * math.cos(theta), radial * math.sin(theta), z)


def frequency_squared(gradient_coefficient, mass_squared, k):
    return gradient_coefficient * math.fsum(x * x for x in k) + mass_squared


def free_energy_density(normalization, gradient_coefficient, mass_squared, temporal, gradient, value):
    return normalization * (temporal * temporal + gradient_coefficient * math.fsum(x * x for x in gradient)
                            + mass_squared * value * value) / 2


def omega(a, mass, k):
    return math.sqrt(frequency_squared(a * a, mass * mass, k))


def boost_mode(energy, k, velocity):
    beta2 = math.fsum(x * x for x in velocity)
    gamma = 1 / math.sqrt(1 - beta2)
    dot = math.fsum(x * y for x, y in zip(k, velocity))
    # gamma²/(gamma+1) equals (gamma-1)/beta² without a small-beta division.
    coefficient = gamma * gamma * dot / (gamma + 1) - gamma * energy
    return gamma * (energy - dot), tuple(x + coefficient * v for x, v in zip(k, velocity))


def principal_coefficients(a, v):
    gamma2 = 1 / (1 - v * v)
    return gamma2 * (1 - a * a * v * v), gamma2 * v * (a * a - 1), gamma2 * (v * v - a * a)


def simpson(function, lower, upper, n):
    step = (upper - lower) / n
    return step * math.fsum((1 if i in (0, n) else 4 if i % 2 else 2) * function(lower + step * i)
                            for i in range(n + 1)) / 3


def bump_derivative(z):
    return -6 * z * (1 - z * z) ** 2 if abs(z) < 1 else 0.0


def laboratory_times(a, v, length=1.0):
    parallel = 2 * a * length * (1 - v * v) / (a * a - v * v)
    transverse = 2 * length * math.sqrt((1 - v * v) / (a * a - v * v))
    return parallel, transverse


def laboratory_events(a, v, length=1.0):
    gamma = 1 / math.sqrt(1 - v * v)
    # Longitudinal lab rod has preferred simultaneous separation L/gamma.
    separation = length / gamma
    forward = separation / (a - v)
    detector_x = v * forward + separation
    # Start leftward ray at the detector; intercept source x=v*t.
    backward = (detector_x - v * forward) / (a + v)
    parallel = (forward + backward) / gamma
    # Transverse: preferred detector/source x=v*t, separation in y remains L.
    first = length / math.sqrt(a * a - v * v)
    second = first
    transverse = (first + second) / gamma
    return parallel, transverse, gamma * forward * (1 - a * v)


def run_checks():
    failures = []
    rng = random.Random(20261004)
    max_derivative_error = max_dispersion_error = max_boost_norm_error = 0.0
    min_boosted_frequency = min_average_mode_energy = math.inf
    superluminal_group_cases = 0
    for i in range(2000):
        a = 1 + 10 ** rng.uniform(-3, 1)
        mass = 0.0 if i % 5 == 0 else 10 ** rng.uniform(-2, 1)
        norm = 10 ** rng.uniform(-3, 3)
        k = tuple(norm * x for x in direction(rng))
        energy = omega(a, mass, k)
        norm = math.sqrt(math.fsum(x * x for x in k))
        phase, group = energy / norm, a * a * norm / energy
        dispersion = a * a * norm * norm + mass * mass
        error = abs(energy * energy - dispersion) / max(1, dispersion)
        max_dispersion_error = max(max_dispersion_error, error)
        average_mode_energy = free_energy_density(1, a * a, mass * mass, energy, k, 1) / 2
        min_average_mode_energy = min(min_average_mode_energy, average_mode_energy)
        if (error > 1e-12 or energy <= 0 or average_mode_energy <= 0
                or not math.isclose(average_mode_energy, energy * energy / 2, rel_tol=1e-12)
                or phase < a * (1 - 1e-12) or group > a * (1 + 1e-12)
                or not math.isclose(phase * group, a * a, rel_tol=1e-12)):
            failures.append(f"Positive-mode dispersion/velocity case {i}")
        expected_superluminal = a * a * (a * a - 1) * norm * norm > mass * mass
        if (group > 1) != expected_superluminal:
            failures.append(f"Group-speed crossover {i}")
        superluminal_group_cases += int(group > 1)
        step = 1e-5 * max(norm, mass / a, 1e-12)
        for j in range(3):
            plus, minus = list(k), list(k)
            plus[j] += step
            minus[j] -= step
            measured = (omega(a, mass, plus) - omega(a, mass, minus)) / (2 * step)
            expected = a * a * k[j] / energy
            max_derivative_error = max(max_derivative_error, abs(measured - expected))
            if not math.isclose(measured, expected, rel_tol=1e-7, abs_tol=1e-7):
                failures.append(f"Dispersion derivative {i},{j}")
        speed = rng.uniform(0, .999)
        velocity = tuple(speed * x for x in direction(rng))
        boosted, boosted_k = boost_mode(energy, k, velocity)
        min_boosted_frequency = min(min_boosted_frequency, boosted)
        measured_norm = boosted * boosted - math.fsum(x * x for x in boosted_k)
        expected_norm = (a * a - 1) * norm * norm + mass * mass
        max_boost_norm_error = max(max_boost_norm_error, abs(measured_norm - expected_norm) / max(1, expected_norm))
        if boosted <= 0 or not math.isclose(measured_norm, expected_norm, rel_tol=1e-9, abs_tol=1e-9):
            failures.append(f"Boosted on-shell mode {i}")

    boost_rows = []
    for a in map(Fraction, ("3/2", "2", "3")):
        for v in map(Fraction, ("0", "1/4", "1/3", "1/2", "2/3", "3/4", "9/10")):
            ptt, ptx, pxx = principal_coefficients(a, v)
            if ptt * pxx - ptx * ptx != -a * a or (ptt > 0) != (v * a < 1):
                failures.append("Exact boosted principal matrix")
            boost_rows.append({"a": str(a), "v": str(v), "Ptt": str(ptt), "Ptx": str(ptx), "Pxx": str(pxx),
                               "slice": "Cauchy" if ptt > 0 else "characteristic" if ptt == 0 else "non_Cauchy"})
    bad_a, bad_v = Fraction(2), Fraction(3, 4)
    bad_ptt, _, _ = principal_coefficients(bad_a, bad_v)
    bad_kx, bad_kperp, bad_mass = Fraction(0), Fraction(1), Fraction(0)
    negative_discriminant = bad_a * bad_a * bad_kx * bad_kx + bad_ptt * (bad_a * bad_a * bad_kperp * bad_kperp + bad_mass * bad_mass)
    if negative_discriminant != Fraction(-80, 7):
        failures.append("Transverse bad-slicing control")

    exact_integral = 72 * (Fraction(1, 3) - Fraction(4, 5) + Fraction(6, 7) - Fraction(4, 9) + Fraction(1, 11))
    if exact_integral != Fraction(1024, 385):
        failures.append("Compact-wave energy integral")
    wave_rows = []
    max_wave_error = max_wave_resolution_difference = max_flux_error = 0.0
    for a in (1.1, 2.0, 10.0):
        for time in (0.0, .3, 1.0, 3.0):
            expected = a * a * float(exact_integral)

            def density(x):
                derivative = bump_derivative(x - a * time)
                return a * a * derivative * derivative

            low, high = (simpson(density, a * time - 1, a * time + 1, n) for n in (256, 512))
            error = max(abs(low - expected), abs(high - expected)) / expected
            resolution = abs(high - low) / expected
            max_wave_error = max(max_wave_error, error)
            max_wave_resolution_difference = max(max_wave_resolution_difference, resolution)
            if error > 1e-8 or resolution > 1e-8:
                failures.append(f"Compact-wave energy {a},{time}")
            for z in (-.75, -.25, .25, .75):
                spatial = bump_derivative(z)
                temporal = -a * spatial
                e = (temporal * temporal + a * a * spatial * spatial) / 2
                flux = -a * a * temporal * spatial
                max_flux_error = max(max_flux_error, abs(flux / e - a))
                if not math.isclose(flux / e, a, rel_tol=1e-12):
                    failures.append("Compact-wave energy flux")
            wave_rows.append({"a": a, "time": time, "energy_exact": expected,
                              "energy_256": low, "energy_512": high})

    max_relay_error = 0.0
    min_relay_ratio = math.inf
    for i in range(2000):
        a = 1 + 10 ** rng.uniform(-3, 1)
        v = rng.uniform(0, .999)
        time = 10 ** rng.uniform(-3, 2)
        delay = 0.0 if i % 5 == 0 else 10 ** rng.uniform(-3, 2)
        gamma = 1 / math.sqrt(1 - v * v)
        received = a * time / (a - v)
        departure = received + gamma * delay * time
        departure_x = v * departure
        arrived = departure + departure_x / a
        expected = (a + v) / (a - v) + gamma * delay * (1 + v / a)
        measured = arrived / time
        max_relay_error = max(max_relay_error, abs(measured - expected) / expected)
        min_relay_ratio = min(min_relay_ratio, measured)
        if measured < 1 or not math.isclose(measured, expected, rel_tol=1e-10, abs_tol=1e-10):
            failures.append(f"Preferred-frame relay {i}")
    if Fraction(4) * (1 - Fraction(9, 10) ** 2) / (2 - Fraction(9, 10)) ** 2 != Fraction(76, 121) or (2 + Fraction(9, 10)) / (2 - Fraction(9, 10)) != Fraction(29, 11):
        failures.append("Exact reciprocal/preferred rule comparison")

    lab_rows = []
    max_lab_error = max_small_v_error = 0.0
    for a in (1.001, 1.1, 2.0, 10.0):
        for v in (0.0, .1, .5, .9, .99):
            closed = laboratory_times(a, v)
            parallel, transverse, forward_tprime = laboratory_events(a, v)
            max_lab_error = max(max_lab_error, *(abs(x - y) / y for x, y in zip((parallel, transverse), closed)))
            if not all(x > 0 and math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-10) for x, y in zip((parallel, transverse), closed)):
                failures.append(f"Laboratory timing {a},{v}")
            lab_rows.append({"a": a, "v": v, "parallel_time_over_L": parallel,
                             "transverse_time_over_L": transverse,
                             "parallel_forward_coordinate_time_over_L": forward_tprime,
                             "tprime_slice_Cauchy": a * v < 1})
        v = .001
        parallel, transverse = laboratory_times(a, v)
        prediction = -(1 - 1 / (a * a)) * v * v / a
        relative = abs((parallel - transverse) - prediction) / abs(prediction)
        max_small_v_error = max(max_small_v_error, relative)
        if relative > 1e-4:
            failures.append(f"Small-v timing expansion {a}")
    for v in (0.0, .1, .5, .9, .99):
        if not all(math.isclose(value, 2, rel_tol=1e-12) for value in laboratory_times(1, v)):
            failures.append(f"Luminal lab timing {v}")

    negative_controls = []
    for description, gradient_coefficient, mass_squared, k, classification in (
            ("m_squared=-1, k=0", 4, -1, (0, 0, 0), "low_k_instability"),
            ("a_squared=-1, m=0, |k|=1", -1, 0, (1, 0, 0), "gradient_instability")):
        squared = frequency_squared(gradient_coefficient, mass_squared, k)
        negative_controls.append({"control": description, "omega_squared": squared,
                                  "growth_rate": math.sqrt(-squared) if squared < 0 else None,
                                  "classification": classification if squared < 0 else "control_not_detected"})
    sign_energies = [free_energy_density(-1, 4, 1, amplitude, (0, 0, 0), amplitude)
                     for amplitude in (1, 2)]
    negative_controls.append({"control": "overall_action_sign=-1", "energy_at_unit_amplitude": sign_energies[0],
                              "energy_at_double_amplitude": sign_energies[1], "classification": "negative_energy_unbounded_scaling"})
    positive_control_energies = [free_energy_density(1, 4, 1, amplitude, (0, 0, 0), amplitude)
                                 for amplitude in (1, 2)]
    controls_detected = (all(row["omega_squared"] < 0 for row in negative_controls[:2])
                         and sign_energies[1] == 4 * sign_energies[0] < sign_energies[0] < 0
                         and all(x > 0 for x in positive_control_energies)
                         and frequency_squared(4, 1, (0, 0, 0)) > 0
                         and frequency_squared(1, 0, (1, 0, 0)) > 0)
    if not controls_detected:
        failures.append("Pathological negative controls were accepted")

    return {"run_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "claim_status": "conditional fixed-background free model; no novelty, viable FTL channel or observational claim",
            "inputs": "deterministic generated modes/worldlines; no private files or observations",
            "parameters": {"seed": 20261004, "a": "1+10**Uniform(-3,1)", "mass": "0 every fifth case; otherwise 10**Uniform(-2,1)",
                           "mode_k_norm": "10**Uniform(-3,3)", "boost_range": [0, .999], "wave_resolutions": [256, 512],
                           "finite_difference_step": "1e-5*max(|k|,m/a,1e-12), registration amendment1"},
            "tolerances": {"dispersion_and_average_energy_relative": 1e-12, "finite_difference_relative_absolute": 1e-7, "boost_norm_relative_absolute": 1e-9,
                           "wave_energy_relative": 1e-8, "worldline_and_lab_relative_absolute": 1e-10,
                           "small_v_relative": 1e-4},
            "passed": not failures, "failures": failures,
            "checks": {"mode_cases": 2000, "derivative_components": 6000, "boosted_modes": 2000,
                       "superluminal_group_modes": superluminal_group_cases, "exact_principal_matrices": len(boost_rows),
                       "compact_wave_cells": len(wave_rows), "preferred_relay_cases": 2000,
                       "laboratory_cells": len(lab_rows), "small_v_cases": 4, "luminal_laboratory_cases": 5,
                       "pathological_controls_detected": controls_detected},
            "residuals": {"max_dispersion_scaled_error": max_dispersion_error, "max_derivative_absolute_error": max_derivative_error,
                          "max_boost_norm_scaled_error": max_boost_norm_error, "min_boosted_frequency": min_boosted_frequency,
                          "min_average_preferred_mode_energy": min_average_mode_energy,
                          "max_wave_energy_relative_error": max_wave_error, "max_wave_resolution_difference": max_wave_resolution_difference,
                          "max_flux_ratio_error": max_flux_error, "max_relay_relative_error": max_relay_error,
                          "min_preferred_relay_ratio": min_relay_ratio, "max_laboratory_relative_error": max_lab_error,
                          "max_small_v_expansion_relative_error": max_small_v_error},
            "principal_matrices": boost_rows, "bad_slicing_transverse_discriminant_over4": str(negative_discriminant),
            "compact_wave_energy": wave_rows, "laboratory_timing": lab_rows, "negative_controls": negative_controls}


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
