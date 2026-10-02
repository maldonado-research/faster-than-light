#!/usr/bin/env python3
"""Internal, independent exact-rational audit of a fixed-u free scalar toy model.

No observational input or parent-agent implementation is imported. Python >=3.10,
standard library only. All checks below are conditional on the stated action.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import random
import sys


def encoded(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encoded(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encoded(v) for v in value]
    return value


def coefficients(a, v, gamma):
    return (gamma**2 * (1-a*a*v*v),
            gamma**2 * v * (a*a-1),
            gamma**2 * (v*v-a*a))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    rng = random.Random(2026100104)
    counts = {}
    failures = []

    def check(kind, condition, detail=''):
        counts[kind] = counts.get(kind, 0) + 1
        if not condition:
            failures.append({'kind': kind, 'detail': detail})

    # Exact rational Lorentz boosts from Pythagorean parametrization.
    for _ in range(400):
        p, q = rng.randint(1, 17), rng.randint(1, 17)
        v = F(p*p-q*q, p*p+q*q)
        gamma = F(p*p+q*q, 2*p*q)
        a = F(rng.randint(11, 55), 10)
        A, B, C = coefficients(a, v, gamma)
        check('lorentz_identity', gamma*gamma*(1-v*v) == 1)
        check('boost_principal_determinant', A*C-B*B == -a*a)
        check('cauchy_sign', (A > 0) == (abs(v)*a < 1))
        check('characteristic_boundary', (A == 0) == (abs(v)*a == 1))

        # Different routes to the transformed symbol: transformed tensor vs
        # directly transformed derivatives of an arbitrary exact covector.
        w = F(rng.randint(-40, 40), 7)
        kx = F(rng.randint(-20, 20), 5)
        ky = F(rng.randint(-20, 20), 5)
        kz = F(rng.randint(-20, 20), 5)
        wp = gamma*(w-v*kx)
        kxp = gamma*(kx-v*w)
        original = w*w-a*a*(kx*kx+ky*ky+kz*kz)
        transformed = A*wp*wp-2*B*wp*kxp+C*kxp*kxp-a*a*(ky*ky+kz*kz)
        check('principal_symbol_invariance', original == transformed)

        # Exact positive-frequency on-shell states: choose frequency first,
        # obtain a nonnegative mass squared, then verify transformed EOM.
        spatial_norm2 = kx*kx+ky*ky+kz*kz
        w = a*(abs(kx)+abs(ky)+abs(kz)+1)
        m2 = w*w-a*a*spatial_norm2
        wp = gamma*(w-v*kx)
        kxp = gamma*(kx-v*w)
        check('preferred_positive_mass_square', m2 >= 0)
        check('boosted_positive_mode_frequency', wp > 0)
        check('minkowski_timelike_on_shell', w*w-spatial_norm2 > 0)
        check('boosted_on_shell_equation',
              A*wp*wp-2*B*wp*kxp+C*kxp*kxp-a*a*(ky*ky+kz*kz) == m2)
        check('boosted_symplectic_flux_identity',
              A*wp-B*kxp == gamma*(w-a*a*v*kx))
        discriminant = a*a*kxp*kxp+A*(a*a*(ky*ky+kz*kz)+m2)
        check('on_shell_discriminant_square', discriminant == (A*wp-B*kxp)**2)

        # Classical translation charge on preferred Cauchy slices. Physical
        # momentum density is -pi*phi_x. This pointwise square completion
        # holds even when t'=constant is not a Cauchy surface.
        pi = F(rng.randint(-20, 20), 3)
        fx = F(rng.randint(-20, 20), 3)
        fy = F(rng.randint(-20, 20), 3)
        fz = F(rng.randint(-20, 20), 3)
        phi = F(rng.randint(-20, 20), 3)
        e = (pi*pi+a*a*(fx*fx+fy*fy+fz*fz)+m2*phi*phi)/2
        q = e+v*pi*fx
        squares = ((pi+v*fx)**2+(a*a-v*v)*fx*fx
                   +a*a*(fy*fy+fz*fz)+m2*phi*phi)/2
        check('translation_charge_square_completion', q == squares)
        check('translation_charge_nonnegative', q >= 0)

        # Moving-boundary finite-propagation estimate: S_x=-a²*pi*phi_x.
        flux = -a*a*pi*fx
        check('finite_propagation_flux_bound', abs(flux) <= a*e)

        # Independent event intersections for moving laboratory round trip.
        # Rest separation L'=1, preferred ray speeds +a and -a, zero delay.
        L = F(1)
        t1 = L/(gamma*(a-v))
        x1 = a*t1
        t2 = 2*a*t1/(a+v)
        x2 = v*t2
        t1p = gamma*(t1-v*x1)
        t2p = gamma*(t2-v*x2)
        outgoing = L*(1-a*v)/(a-v)
        incoming = L*(1+a*v)/(a+v)
        check('round_trip_ray_intersections',
              x1 == v*t1+L/gamma and x1-a*(t2-t1) == x2)
        check('round_trip_preferred_future', t1 > 0 and t2 > t1)
        check('one_way_lab_formula', t1p == outgoing)
        check('round_trip_lab_formula', t2p == outgoing+incoming
              == 2*a*L*(1-v*v)/(a*a-v*v))
        check('round_trip_lab_future', t2p > 0)
        # Transverse baseline L'=1 is unchanged by an x boost. Intersections
        # obey (a²-v²)t1²=1, identical for the return leg. Squared proper
        # round-trip time avoids introducing irrational floating values.
        transverse_preferred_t1_squared = L*L/(a*a-v*v)
        transverse_proper_t2_squared = 4*transverse_preferred_t1_squared/(gamma*gamma)
        check('transverse_direct_intersection',
              a*a*transverse_preferred_t1_squared
              == v*v*transverse_preferred_t1_squared+L*L)
        check('transverse_proper_round_trip',
              transverse_proper_t2_squared == 4*L*L*(1-v*v)/(a*a-v*v))
        squared_orientation_ratio = t2p*t2p/transverse_proper_t2_squared
        check('orientation_ratio',
              squared_orientation_ratio == a*a*(1-v*v)/(a*a-v*v))
        check('orientation_sign',
              (squared_orientation_ratio < 1) == (v != 0))

    # Reproducible counterexample to treating a negative boosted kinetic
    # coefficient as evidence of a physical ghost.
    a, v, gamma = F(2), F(4, 5), F(5, 3)
    A, B, C = coefficients(a, v, gamma)
    # Arbitrary t' spatial mode kx'=0, kperp'=1, m=0 has imaginary frequency.
    disc = A*a*a
    imaginary_frequency_square = a*a/A
    check('noncauchy_transverse_illposedness', disc < 0
          and imaginary_frequency_square == -F(12, 13))
    # The real preferred on-shell wave (omega,kx)=(2,1) transforms to (2,-1),
    # with positive observer frequency despite A<0.
    wp = gamma*(F(2)-v)
    kxp = gamma*(F(1)-v*F(2))
    check('negative_coefficient_positive_mode_counterexample',
          A == -F(13, 3) and wp == 2 and kxp == -1
          and A*wp*wp-2*B*wp*kxp+C*kxp*kxp == 0)
    check('positive_frequency_negative_noncauchy_flux',
          A*wp-B*kxp == -2 and wp > 0)
    # At kx'=-1 both roots represent preferred positive-frequency waves.
    roots = [F(14, 13), F(2)]
    check('folded_spatial_spectrum_roots',
          all(A*r*r+2*B*r+C == 0 and gamma*(r-v) > 0 for r in roots))
    # The same example has an Einstein-synchronized outbound reception
    # time before emission, but a positive emitter's proper round-trip time.
    out = (1-a*v)/(a-v)
    back = (1+a*v)/(a+v)
    check('negative_one_way_positive_round_trip',
          out == -F(1, 2) and back == F(13, 14) and out+back == F(3, 7))

    # Source-free wave packets are stable in preferred time even where group
    # speed is >1. An exact massive example avoids radicals: a=2,k=3,m²=28,
    # omega=8, phase=8/3, group=3/2, front=2, covector timelike.
    check('massive_superluminal_group_example',
          F(8)**2 == F(2)**2*F(3)**2+28 and F(4)*3/8 == F(3, 2))

    result = {
        'status': 'PASS' if not failures else 'FAIL',
        'run_utc': datetime.now(timezone.utc).isoformat(),
        'evidence_tier': 'internal independent exact-rational synthetic checks',
        'action_assumptions': {
            'background': 'Minkowski eta=(-,+,+,+), fixed constant unit timelike u',
            'preferred_equation': 'phi_tt-a^2 Laplacian(phi)+m^2 phi=J',
            'a': '>1', 'm_squared': '>=0', 'c': 1,
            'boundary': 'free field or retarded source; decay/finite energy for charges',
            'physical_couplings': 'unspecified; no empirical FTL claim',
        },
        'seed': 2026100104,
        'python': platform.python_version(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'counts': counts,
        'total_checks': sum(counts.values()),
        'failures': failures,
        'exact_noncauchy_example': {
            'a': a, 'v': v, 'gamma': gamma, 'A': A, 'B': B, 'C': C,
            'transverse_discriminant': disc,
            'arbitrary_transverse_omega_prime_squared': imaginary_frequency_square,
            'preferred_positive_mode': {'omega': F(2), 'kx': F(1)},
            'boosted_same_mode': {'omega_prime': wp, 'kx_prime': kxp},
            'two_positive_roots_at_kx_prime_minus_one': roots,
            'one_way_t_prime': out, 'return_delta_t_prime': back,
            'round_trip_emitter_proper_time': out+back,
        },
        'limitations': [
            'No observational inputs, production mechanism, detector model, or interacting completion.',
            'Finite retarded support and no-loop statement additionally have analytic proofs in AUDIT.md.',
            'Coordinate covariance transforms u; a fixed-u state is not invariant under boosts.',
            'Checks are internal assistance, not external peer review or physical replication.',
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(encoded(result), indent=2)+'\n')
    print(f"{result['status']}: {result['total_checks']} exact checks; output={args.out}")
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
