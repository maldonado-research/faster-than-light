#!/usr/bin/env python3
"""Registered 1+1D coupled-scalar retarded-kernel checks; stdlib/Python >=3.10.

Deterministic mathematical parameters, no observations or physical detector model.
Run after reading REGISTRATION.md: python3 reproduce.py --out RESULTS.json
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import platform


TERMS = 12


def j0_squared(z2: float) -> float:
    if not 0 <= z2 <= 1:
        raise ValueError("Registered Bessel argument square must be in [0,1]")
    terms = [1.0]
    term = 1.0
    for j in range(1, TERMS):
        term *= -z2 / (4*j*j)
        terms.append(term)
    return math.fsum(terms)


def two_j1_over_z_squared(z2: float) -> float:
    if not 0 <= z2 <= 1:
        raise ValueError("Registered Bessel argument square must be in [0,1]")
    terms = [1.0]
    term = 1.0
    for j in range(1, TERMS):
        term *= -z2 / (4*j*(j+1))
        terms.append(term)
    return math.fsum(terms)


def weight(i: int, n: int) -> int:
    return 1 if i in (0, n) else 4 if i % 2 else 2


def triangle(a: float, m: float, big_m: float, t: float, r: float, n: int):
    """Return two normalized integrals on the exact finite intercone triangle."""
    if not (a > 1 and t < r < a*t and n > 0 and n % 2 == 0):
        raise ValueError("Invalid registered triangle domain or Simpson resolution")
    gap = a*t-r
    delta = a*a-1
    alpha = (a+1)/(a-1)
    beta = (a-1)/(a+1)
    rows_h, rows_k = [], []
    maximum_argument_square = 0.0
    for i in range(n+1):
        u = i/n
        remaining = 1-u
        row_h, row_k = [], []
        for j in range(n+1):
            z = j/n
            v = remaining*z
            # Compute 1-u-v without cancellation at the triangle boundary.
            w = remaining*(1-z)
            fast_square = gap/a*w*(2*t-gap/a*(1+alpha*u+beta*v))
            slow_square = 4*gap*gap*u*v/delta
            fast_z2, slow_z2 = m*m*fast_square, big_m*big_m*slow_square
            maximum_argument_square = max(maximum_argument_square, fast_z2, slow_z2)
            fast_j0 = j0_squared(fast_z2)
            slow_j0 = j0_squared(slow_z2)
            slow_factor = two_j1_over_z_squared(slow_z2)
            jacobian = remaining
            wj = weight(j, n)
            row_h.append(wj*jacobian*fast_j0*slow_j0)
            row_k.append(wj*jacobian*u*v*fast_j0*slow_factor)
        wi = weight(i, n)
        rows_h.append(wi*math.fsum(row_h))
        rows_k.append(wi*math.fsum(row_k))
    denominator = 9*n*n
    return (math.fsum(rows_h)/denominator,
            math.fsum(rows_k)/denominator,
            maximum_argument_square)


def moment_polynomial(p: int, q: int) -> Q:
    # Integrate v first, then expand (1-u)^(q+1) and integrate u exactly.
    return sum((Q((-1)**j*math.comb(q+1, j), (q+1)*(p+j+1))
                for j in range(q+2)), Q(0))


def exact_front_data(a: Q, m: Q, big_m: Q, g: Q, t: Q, r: Q):
    gap = a*t-r
    delta = a*a-1
    h0 = gap*gap/(4*a*delta)
    k0 = gap**4/(48*a*delta*delta)
    error_h = m*m*t*gap/(2*a)+big_m*big_m*gap*gap/(4*delta)
    error_k = m*m*t*gap/(2*a)+big_m*big_m*gap*gap/(8*delta)
    small_time_factor = (1-m*m*t*t/4)*(1-big_m*big_m*t*t/4)**2
    tail_first = g**4*t*gap**7/(20160*a*a*delta**3)
    tail_ratio = g*g*t*gap**3/(2880*a*delta)
    tail_upper = tail_first/(1-tail_ratio)
    return {
        "gap": gap, "H0": h0, "K0": k0,
        "relative_error_bound_H": error_h, "relative_error_bound_K": error_k,
        "small_time_factor": small_time_factor,
        "chi_exact_lower": g*g*small_time_factor*k0,
        "chi_exact_refined_lower": g*g*(1-error_k)*k0,
        "chi_exact_upper": g*g*k0+tail_upper,
        "higher_mixing_upper": tail_upper, "higher_mixing_geometric_ratio": tail_ratio,
    }


def run_checks():
    failures = []
    counts = {"exact_triangle_moments": 0, "transform_cases": 0,
              "exact_transform_identities": 0, "massive_front_cell_evaluations": 0,
              "massive_front_quadrature_grids": 0, "massive_front_scalar_integrals": 0,
              "massless_quadrature_grids": 0, "massless_normalization_checks": 0,
              "full_response_bound_cells": 0, "front_monotonicity_checks": 0,
              "front_power_checks": 0, "front_slope_checks": 0,
              "stability_controls": 0}

    def require(condition, label):
        if not condition:
            failures.append(label)

    for p in range(6):
        for q in range(6):
            actual = moment_polynomial(p, q)
            expected = Q(math.factorial(p)*math.factorial(q), math.factorial(p+q+2))
            counts["exact_triangle_moments"] += 1
            require(actual == expected, f"Triangle moment {p},{q}")

    for a in (Q(5,4), Q(3,2), Q(2), Q(3)):
        for laplace in (Q(1,2), Q(1), Q(2)):
            for k in (Q(0), Q(1)):
                delta = a*a-1
                fast = laplace*laplace+a*a*k*k
                slow = laplace*laplace+k*k
                h_original = 1/(fast*slow)
                h_partial = (a*a/fast-1/slow)/(delta*laplace*laplace)
                k_original = 1/(fast*slow*slow)
                k_partial = (a**4/(delta**2*laplace**4*fast)
                             -a*a/(delta**2*laplace**4*slow)
                             -1/(delta*laplace*laplace*slow*slow))
                counts["transform_cases"] += 1
                counts["exact_transform_identities"] += 2
                require(h_original == h_partial, "Exact H transform identity")
                require(k_original == k_partial, "Exact K transform identity")

    massless_rows = []
    for n in (32, 64):
        h, k, argument = triangle(2., 0., 0., 1., 1.5, n)
        counts["massless_quadrature_grids"] += 1
        counts["massless_normalization_checks"] += 2
        require(abs(h-.5) <= 1e-13, "Massless H triangle normalization")
        require(abs(k-1/24) <= 1e-13, "Massless K triangle normalization")
        massless_rows.append({"resolution": n, "H_integral": h, "K_integral": k,
                              "maximum_bessel_argument_square": argument})

    a, m, big_m, g = Q(2), Q(1), Q(1), Q(1,2)
    selected = []
    for t in (Q(1,100), Q(1,20), Q(1,10)):
        for ratio in (Q(5,4), Q(3,2), Q(7,4)):
            selected.append(("stable", t, ratio*t))
    for gap_ratio in (Q(1,2), Q(1,4), Q(1,8), Q(1,16), Q(1,32)):
        t = Q(1,10)
        selected.append(("front", t, a*t-gap_ratio*t))

    rows = []
    maximum_resolution_change = 0.0
    for kind, t, r in selected:
        exact = exact_front_data(a, m, big_m, g, t, r)
        require(g*g <= m*m*big_m*big_m and m*t <= 1 and big_m*t <= 1,
                "Stable small-time parameters")
        require(t < r < a*t, "Intercone parameter")
        require(exact["chi_exact_lower"] > 0, "Strictly positive full chi lower bound")
        require(exact["chi_exact_refined_lower"] > 0, "Refined chi lower bound")
        require(exact["higher_mixing_geometric_ratio"] < 1, "Volterra geometric bound")
        require(exact["chi_exact_lower"] <= exact["chi_exact_upper"], "Full chi interval")
        counts["full_response_bound_cells"] += 1
        values = []
        for n in (64, 128):
            hi, ki, max_argument = triangle(float(a), float(m), float(big_m),
                                             float(t), float(r), n)
            values.append((hi, ki, max_argument))
            counts["massive_front_quadrature_grids"] += 1
            counts["massive_front_scalar_integrals"] += 2
        h_low, k_low, _ = values[0]
        h_high, k_high, maximum_argument = values[1]
        h_change, k_change = abs(h_high-h_low)/h_high, abs(k_high-k_low)/k_high
        maximum_resolution_change = max(maximum_resolution_change, h_change, k_change)
        require(h_change <= 1e-9 and k_change <= 1e-9, "Two-resolution agreement")
        h_ratio, k_ratio = 2*h_high, 24*k_high
        require(1-float(exact["relative_error_bound_H"])-2e-12 <= h_ratio <= 1+2e-12,
                "Analytic massive H ratio bound")
        require(1-float(exact["relative_error_bound_K"])-2e-12 <= k_ratio <= 1+2e-12,
                "Analytic massive K ratio bound")
        counts["massive_front_cell_evaluations"] += 1
        rows.append({
            "case": kind, "t": str(t), "r": str(r), "r_over_t": str(r/t),
            "gap": str(exact["gap"]), "resolutions": [64,128],
            "H_integrals": [h_low,h_high], "K_integrals": [k_low,k_high],
            "relative_resolution_change_H": h_change, "relative_resolution_change_K": k_change,
            "H_over_H0": h_ratio, "K_over_K0": k_ratio,
            "H_massive_coefficient": float(exact["H0"])*h_ratio,
            "K_massive_coefficient": float(exact["K0"])*k_ratio,
            "leading_chi_coefficient_g2_K": float(g*g*exact["K0"])*k_ratio,
            "maximum_bessel_argument_square": maximum_argument,
            "analytic_exact_rational_bounds": {name: str(value) for name,value in exact.items()},
            "analytic_bounds_float": {name: float(value) for name,value in exact.items()},
        })

    front = [row for row in rows if row["case"] == "front"]
    for left, right in zip(front, front[1:]):
        for key in ("H_over_H0", "K_over_K0"):
            counts["front_monotonicity_checks"] += 1
            require(left[key] < right[key] < 1, "Front ratio approaches 1")
        for key, power in (("H_massive_coefficient",2), ("K_massive_coefficient",4)):
            observed = left[key]/right[key]
            counts["front_power_checks"] += 1
            require(abs(observed/(2**power)-1) <= 1e-2, "Near-front power suppression")
    final = front[-1]
    gap = float(Q(final["gap"]))
    slope_rows = []
    for key, expected in (("H_over_H0",float(m*m*Q(1,10)/(6*a))),
                          ("K_over_K0",float(m*m*Q(1,10)/(10*a)))):
        actual = (1-final[key])/gap
        relative = abs(actual-expected)/expected
        counts["front_slope_checks"] += 1
        require(relative <= 1e-2, "Registered leading mass slope")
        slope_rows.append({"ratio":key,"actual_slope":actual,
                           "expected_slope":expected,"relative_error":relative})

    controls = []
    for mass, mixing, expected_stable in ((Q(1),Q(1,2),True),
                                          (Q(0),Q(1,2),False),
                                          (Q(0),Q(-1,2),False)):
        # At k=0 and m=M the exact eigenvalues are mass^2 +/- |g|.
        minimum = mass*mass-abs(mixing)
        determinant = mass**4-mixing*mixing
        is_stable = determinant >= 0 and minimum >= 0
        counts["stability_controls"] += 1
        require(is_stable == expected_stable, "Stable/full-massless control")
        controls.append({"m_and_M":str(mass),"g":str(mixing),
                         "minimum_omega2_at_k0":str(minimum),
                         "potential_determinant":str(determinant),"stable":is_stable})

    expected_counts = {"exact_triangle_moments":36,"transform_cases":24,
                       "exact_transform_identities":48,"massive_front_cell_evaluations":14,
                       "massive_front_quadrature_grids":28,"massive_front_scalar_integrals":56,
                       "massless_quadrature_grids":2,"massless_normalization_checks":4,
                       "full_response_bound_cells":14,"front_monotonicity_checks":8,
                       "front_power_checks":8,"front_slope_checks":2,"stability_controls":3}
    require(counts == expected_counts, "Registered counts not completed")
    folder = Path(__file__).resolve().parent
    return {
        "run_utc":datetime.now(timezone.utc).isoformat(),"python":platform.python_version(),
        "dependencies":"Python >=3.10 standard library only",
        "claim_status":"Conditional 1+1D scalar analogue; no observations, detector prediction or novelty claim",
        "prior_formula_knowledge":"Algebraic derivation disclosed before this decisive run",
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "registration_sha256":hashlib.sha256((folder/"REGISTRATION.md").read_bytes()).hexdigest(),
        "root_registration_sha256":hashlib.sha256((folder.parent/"REGISTRATION.md").read_bytes()).hexdigest(),
        "passed":not failures,"failures":failures,"checks":counts,
        "parameters":{"a":"2","m":"1","M":"1","g":"1/2","resolutions":[64,128],
                      "bessel_series_terms":TERMS,"distinct_massive_front_cells":len({(t,r) for _,t,r in selected})},
        "tolerances":{"two_resolution_relative":1e-9,"analytic_normalized_bound_slack":2e-12,
                      "massless_normalization_absolute":1e-13,"front_power_and_slope_relative":1e-2},
        "bessel_truncation_bounds_at_z2_le_1":{
            "J0":str(Q(1,4**TERMS*math.factorial(TERMS)**2)),
            "2J1_over_z":str(Q(1,4**TERMS*math.factorial(TERMS)*math.factorial(TERMS+1)))},
        "maximum_relative_resolution_change":maximum_resolution_change,
        "massless_normalization":massless_rows,"cells":rows,
        "final_front_slope_checks":slope_rows,"stability_controls":controls,
        "full_response_status":"Exact positive lower/upper bounds from convergent Volterra series; full kernel not numerically Fourier-inverted",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    output = args.out.resolve()
    protected = {
        Path(__file__).resolve(), folder / "REGISTRATION.md",
        folder.parent / "REGISTRATION.md", folder / "RESULTS.json",
        folder / "INITIAL_REPRODUCE.py", folder / "INITIAL_RESULTS.json",
        folder / "DERIVATION.md", folder / "AUDIT.md",
    }
    if args.out.is_symlink() or output.exists() or output in protected:
        raise SystemExit("Output must be a fresh file and cannot replace protected source or records")
    try:
        result = run_checks()
    except Exception as error:
        result = {"run_utc":datetime.now(timezone.utc).isoformat(),
                  "passed":False,"failures":[f"{type(error).__name__}: {error}"],
                  "unexpected_exception":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x",encoding="utf-8") as stream:
        stream.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({key:result[key] for key in
                      ("passed","failures","checks","maximum_relative_resolution_change")
                      if key in result}))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
