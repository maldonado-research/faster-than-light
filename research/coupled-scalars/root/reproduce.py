#!/usr/bin/env python3
"""Registered, synthetic coupled-scalar checks; Python 3.10+ standard library.

Require an explicit output and refuse overwrite. The tested continuum model is
an assumed scalar analogue, not an observed channel or source/detector system.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path
import platform
import random
import sys
import traceback

FLOAT_TOL = 1e-10
ASYMPTOTIC_TOL = 1e-4
REQUIRED = {'modes': 2000, 'partial_fractions': 200, 'piecewise_kernels': 120,
            'spatial_integrals': 60, 'fast_residue': 20}
EXTRA = {'mass_boundaries': 16, 'coincident_cones': 24,
         'zero_wavenumber_transforms': 6, 'pathological_controls': 4}


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


class Checks:
    def __init__(self):
        self.sections = {}
        self.failures = []

    def case(self, section, label, inputs, conditions, metrics=None):
        row = self.sections.setdefault(section, {'cases': 0, 'passed': 0,
                                                 'failed': 0, 'maxima': {}})
        row['cases'] += 1
        failed = [name for name, value in conditions.items() if not value]
        for name, value in (metrics or {}).items():
            value = float(value)
            if not math.isfinite(value):
                failed.append(name + ': nonfinite metric')
            row['maxima'][name] = max(row['maxima'].get(name, -math.inf), value)
        row['failed' if failed else 'passed'] += 1
        if failed:
            # Retain all inputs and named equations for every failed case.
            self.failures.append(encode({'section': section, 'label': label,
                                         'inputs': inputs, 'failed_checks': failed,
                                         'metrics': metrics or {}}))


def q(rng):
    return F(rng.randint(-20, 20), rng.randint(1, 8))


def matrix(a2, m, M, g, k2):
    return a2*k2 + m*m, k2 + M*M, g


def hamiltonian(a2, m, M, g, k2, fields, velocities, kinetic=(F(1), F(1))):
    phi, chi = fields
    vphi, vchi = velocities
    return (kinetic[0]*vphi*vphi + kinetic[1]*vchi*vchi +
            a2*k2*phi*phi + k2*chi*chi + m*m*phi*phi + M*M*chi*chi)/2 + g*phi*chi


def eigenvalues(A, B, g):
    # Only the larger root uses addition; the smaller uses the exact determinant.
    trace = float(A+B)
    D = math.hypot(float(A-B), 2*float(g))
    high = (trace + D)/2
    determinant = A*B-g*g
    low = float(determinant)/high if high else 0.0
    return high, low, D


def mode_checks(a, m, M, g, k, fields, velocities, beta):
    a2 = a*a
    k2 = sum((ki*ki for ki in k), F(0))
    A, B, off = matrix(a2, m, M, g, k2)
    determinant = A*B-off*off
    mass_det = m*m*M*M-g*g
    determinant_expanded = a2*k2*k2 + (a2*M*M+m*m)*k2 + mass_det
    energy = hamiltonian(a2, m, M, g, k2, fields, velocities)
    phi, chi = fields
    matrix_energy = (sum((v*v for v in velocities), F(0)) +
                     phi*(A*phi+off*chi) + chi*(off*phi+B*chi))/2
    shifted_det = (A-k2)*(B-k2)-off*off
    shifted_expanded = (a2-1)*M*M*k2 + mass_det
    conditions = {
        'mass determinant nonnegative': mass_det >= 0,
        'full matrix energy equals density': matrix_energy == energy,
        'Hamiltonian nonnegative': energy >= 0,
        'exact determinant expansion': determinant == determinant_expanded,
        'exact trace equals diagonal sum': A+B == (a2+1)*k2+m*m+M*M,
        'shifted diagonal nonnegative': A-k2 >= 0 and B-k2 >= 0,
        'shifted determinant nonnegative': shifted_det >= 0,
        'exact shifted determinant expansion': shifted_det == shifted_expanded,
    }
    high, low, D = eigenvalues(A, B, off)
    theta = math.atan2(2*float(off), float(A-B))/2
    cs, sn = math.cos(theta), math.sin(theta)
    vectors = ((cs, sn), (-sn, cs))
    Af, Bf, gf = float(A), float(B), float(off)
    trace, det = float(A+B), float(determinant)
    char_scale = max(1.0, abs(trace)**2, abs(det))
    eig_raw = eig_normalized = char_raw = char_normalized = 0.0
    speed_violation = 0.0
    beta2 = sum((v*v for v in beta), F(0))
    beta_dot = sum((v*ki for v, ki in zip(beta, k)), F(0))
    conditions['ordinary boost subluminal exactly'] = beta2 < 1
    conditions['exact boost Cauchy-Schwarz bound'] = beta_dot*beta_dot <= beta2*k2
    boosted = []
    for label, value, vector in zip(('high', 'low'), (high, low), vectors):
        u, v = vector
        raw = max(abs(Af*u+gf*v-value*u), abs(gf*u+Bf*v-value*v))
        scale = max(1.0, abs(Af), abs(Bf), abs(gf), abs(value))
        normalized = raw/scale
        characteristic = abs(value*value-trace*value+det)
        eig_raw, eig_normalized = max(eig_raw, raw), max(eig_normalized, normalized)
        char_raw = max(char_raw, characteristic)
        char_normalized = max(char_normalized, characteristic/char_scale)
        conditions[label + ' eigenproblem'] = normalized <= FLOAT_TOL
        conditions[label + ' characteristic polynomial'] = characteristic/char_scale <= FLOAT_TOL
        conditions[label + ' nonnegative frequency squared'] = value >= 0
        violation = max(0.0, float(k2)-value)/max(1.0, abs(value), float(k2))
        speed_violation = max(speed_violation, violation)
        conditions[label + ' omega at least |k|'] = violation <= FLOAT_TOL
        omega = math.sqrt(value) if value >= 0 else math.nan
        omega_prime = (omega-float(beta_dot))/math.sqrt(1-float(beta2))
        boosted.append(omega_prime)
        conditions[label + ' ordinary boosted frequency'] = (
            omega_prime == 0.0 if omega == 0.0 else omega_prime > 0.0)
    orthogonal_raw = max(abs(cs*cs+sn*sn-1), abs(cs*(-sn)+sn*cs))
    reconstruction = (high*cs*cs+low*sn*sn,
                      high*sn*sn+low*cs*cs,
                      (high-low)*cs*sn)
    reconstruct_raw = max(abs(reconstruction[0]-Af), abs(reconstruction[1]-Bf),
                          abs(reconstruction[2]-gf))
    reconstruct_norm = reconstruct_raw/max(1.0, abs(Af), abs(Bf), abs(gf), high)
    trace_raw = abs(high+low-trace)
    trace_norm = trace_raw/max(1.0, abs(trace))
    product_raw = abs(high*low-det)
    product_norm = product_raw/char_scale
    conditions.update({'orthonormal eigenvectors': orthogonal_raw <= FLOAT_TOL,
                       'full matrix reconstruction': reconstruct_norm <= FLOAT_TOL,
                       'eigenvalue trace': trace_norm <= FLOAT_TOL,
                       'eigenvalue product': product_norm <= FLOAT_TOL})
    metrics = {'eigen_raw': eig_raw, 'eigen_normalized': eig_normalized,
               'characteristic_raw': char_raw, 'characteristic_normalized': char_normalized,
               'reconstruction_raw': reconstruct_raw, 'reconstruction_normalized': reconstruct_norm,
               'trace_raw': trace_raw, 'trace_normalized': trace_norm,
               'product_raw': product_raw, 'product_normalized': product_norm,
               'orthogonality': orthogonal_raw, 'speed_bound_violation_normalized': speed_violation}
    return conditions, metrics, {'k2': k2, 'energy': energy, 'mass_det': mass_det,
                                  'lambda_high': high, 'lambda_low': low,
                                  'ordinary_boosted_frequencies': boosted}


def kernels(t, x, a):
    if t <= 0:
        return F(0), F(0)
    r = abs(x)
    qa, q1 = max(t-r/a, F(0)), max(t-r, F(0))
    d = a*a-1
    H = (a*qa*qa-q1*q1)/(4*d)
    K = (2*a**3*qa**4-(3*a*a-1)*q1**4-4*d*r*q1**3)/(96*d*d)
    return H, K


def region_kernels(t, r, a):
    # Independently expanded nonsingular inner and fast-front polynomials.
    if t <= 0 or r > a*t:
        return F(0), F(0)
    if r <= t:
        return ((a*t*t-r*r)/(4*a*(a+1)),
                ((2*a+1)*t**4-6*t*t*r*r+(a+2)*r**4/a)/(96*(a+1)**2))
    d = a*a-1
    return ((a*t-r)**2/(4*a*d), (a*t-r)**4/(48*a*d*d))


def equal_speed_kernels(t, x):
    # Null-coordinate retarded rectangle convolutions, not a singular a=1 substitution.
    u, v = (t+x)/2, (t-x)/2
    if t <= 0 or u < 0 or v < 0:
        return F(0), F(0)
    return u*v/2, u*u*v*v/8


def leading_responses(g, H, K):
    return -g*H, g*g*K


def pscale(coeffs, value):
    return [c*value for c in coeffs]


def padd(*polys):
    out = [F(0)]*max(map(len, polys))
    for poly in polys:
        for index, coefficient in enumerate(poly):
            out[index] += coefficient
    return out


def shifted_power(t, speed, power):
    return [F(math.comb(power, index))*t**(power-index)*(-1/speed)**index
            for index in range(power+1)]


def kernel_polynomials(t, a):
    d = a*a-1
    H_fast = pscale(shifted_power(t, a, 2), a/(4*d))
    H_inner = padd(H_fast, pscale(shifted_power(t, F(1), 2), -1/(4*d)))
    K_fast = pscale(shifted_power(t, a, 4), 2*a**3/(96*d*d))
    K_inner = padd(K_fast,
                   pscale(shifted_power(t, F(1), 4), -(3*a*a-1)/(96*d*d)),
                   [F(0)] + pscale(shifted_power(t, F(1), 3), -4*d/(96*d*d)))
    return H_inner, H_fast, K_inner, K_fast


def integrate(poly, lo, hi):
    return sum((c*(hi**(index+1)-lo**(index+1))/(index+1)
                for index, c in enumerate(poly)), F(0))


def evaluate(poly, x):
    total = F(0)
    for coefficient in reversed(poly):
        total = total*x+coefficient
    return total


def run(checks):
    diagnostics = {'interior_grid_evaluations': 0, 'residue_doubling_ratios': 0,
                   'random_mode_saturation_cases': 0, 'random_mode_zero_mixing_cases': 0,
                   'random_mode_zero_frequency_modes': 0}
    rng = random.Random(20261005)
    for index in range(REQUIRED['modes']):
        a = F(rng.randint(5, 20), 4)
        m, M = (F(rng.randint(1, 12), rng.randint(1, 6)) for _ in range(2))
        p = F(rng.randint(-12, 12), 12)
        g = p*m*M
        k = [q(rng) for _ in range(3)]
        fields, velocities = [q(rng) for _ in range(2)], [q(rng) for _ in range(2)]
        beta = [F(rng.randint(-4, 4), 10) for _ in range(3)]
        conditions, metrics, details = mode_checks(a, m, M, g, k, fields, velocities, beta)
        checks.case('modes', index, locals_inputs(a=a, m=m, M=M, g=g, k=k,
                    fields=fields, velocities=velocities, beta=beta), conditions, metrics)
        diagnostics['random_mode_saturation_cases'] += int(abs(p) == 1)
        diagnostics['random_mode_zero_mixing_cases'] += int(g == 0)
        diagnostics['random_mode_zero_frequency_modes'] += sum(
            int(details[name] == 0) for name in ('lambda_high', 'lambda_low'))

    rng = random.Random(20261006)
    for index in range(REQUIRED['partial_fractions']):
        u, v = (F(rng.randint(1, 20), rng.randint(1, 8)) for _ in range(2))
        if u == v:
            v += F(1, 7)
        s, k = (F(rng.randint(1, 12), rng.randint(1, 7)) for _ in range(2))
        Pa, Pb, delta = s*s+u*u*k*k, s*s+v*v*k*k, (u*u-v*v)*k*k
        original_H, original_K = 1/(Pa*Pb), 1/(Pa*Pb*Pb)
        expanded_H = (1/Pb-1/Pa)/delta
        expanded_K = 1/(delta*Pb*Pb)+(1/Pa-1/Pb)/(delta*delta)
        checks.case('partial_fractions', index, locals_inputs(u=u, v=v, s=s, k=k),
                    {'distinct positive speeds': u > 0 and v > 0 and u != v,
                     'positive transform variables': s > 0 and k > 0,
                     'H exact partial fractions': original_H == expanded_H,
                     'K exact partial fractions': original_K == expanded_K})

    for index, s in enumerate(map(F, ('1/4', '1/2', '1', '2', '3', '5'))):
        Pa, Pb = s*s+F(4)*0, s*s+F(1)*0
        checks.case('zero_wavenumber_transforms', index, {'s': s, 'k': 0},
                    {'H nonsingular limit': 1/(Pa*Pb) == 1/s**4,
                     'K nonsingular limit': 1/(Pa*Pb*Pb) == 1/s**6})

    speeds = list(map(F, ('5/4', '3/2', '2', '3', '5')))
    regions = ('nonpositive', 'outside', 'fast_front', 'intercone', 'ordinary_front', 'inner')
    for index in range(10):
        a, t = speeds[index % len(speeds)], F(index+1, index % 4+1)
        for region in regions:
            ti = -t if index % 2 else F(0)
            r = F(index+1, 3)
            if region != 'nonpositive':
                ti = t
                r = {'outside': a*t+t/3, 'fast_front': a*t,
                     'intercone': (a+1)*t/2, 'ordinary_front': t, 'inner': t/3}[region]
            for sign in (-1, 1):
                x = sign*r
                H, K = kernels(ti, x, a)
                expected_H, expected_K = region_kernels(ti, r, a)
                cross_p, self_p = leading_responses(F(1, 2), H, K)
                cross_m, self_m = leading_responses(F(-1, 2), H, K)
                conditions = {'H region polynomial': H == expected_H,
                              'K region polynomial': K == expected_K,
                              'spatial evenness': (H, K) == kernels(ti, -x, a),
                              'H and K nonnegative': H >= 0 and K >= 0,
                              'odd cross Born coefficient': cross_p == -cross_m,
                              'even self Born coefficient': self_p == self_m}
                if region in ('nonpositive', 'outside', 'fast_front'):
                    conditions['exact zero outside/on fast onset'] = H == K == 0
                if region == 'intercone':
                    d = a*a-1
                    conditions['quadratic intercone onset'] = H == (a*ti-r)**2/(4*a*d) and H > 0
                    conditions['quartic intercone onset'] = K == (a*ti-r)**4/(48*a*d*d) and K > 0
                checks.case('piecewise_kernels', f'{index}:{region}:{sign}',
                            locals_inputs(a=a, t=ti, x=x, region=region), conditions)

    for t in map(F, ('1/10', '1/5', '1/2', '1', '2', '3')):
        for ratio in map(F, ('0', '1/3', '1', '4/3')):
            r = ratio*t
            actual = equal_speed_kernels(t, r)
            independent = region_kernels(t, r, F(1))
            checks.case('coincident_cones', f'{t}:{ratio}', {'t': t, 'r': r, 'a': 1},
                        {'null convolution equals nonsingular polynomial': actual == independent,
                         'evenness': actual == equal_speed_kernels(t, -r),
                         'nonnegative': all(v >= 0 for v in actual),
                         'outside zero': actual == (F(0), F(0)) if r >= t else all(v > 0 for v in actual)})

    rng = random.Random(20261008)
    for index in range(REQUIRED['spatial_integrals']):
        a, t = F(rng.randint(5, 20), 4), F(rng.randint(1, 12), rng.randint(1, 8))
        Hi, Hf, Ki, Kf = kernel_polynomials(t, a)
        integrated_H = 2*(integrate(Hi, F(0), t)+integrate(Hf, t, a*t))
        integrated_K = 2*(integrate(Ki, F(0), t)+integrate(Kf, t, a*t))
        conditions = {'exact H spatial normalization': integrated_H == t**3/6,
                      'exact K spatial normalization': integrated_K == t**5/120}
        for j in range(1, 9):
            for region, r, hp, kp in (
                    ('inner', j*t/9, Hi, Ki),
                    ('intercone', t+j*(a-1)*t/9, Hf, Kf)):
                H, K = kernels(t, r, a)
                conditions[f'{region} grid {j} polynomial agreement'] = (
                    H == evaluate(hp, r) and K == evaluate(kp, r))
                conditions[f'{region} grid {j} nonnegativity'] = H >= 0 and K >= 0
                diagnostics['interior_grid_evaluations'] += 1
        checks.case('spatial_integrals', index, {'a': a, 't': t}, conditions,
                    {'H_exact_residual': abs(integrated_H-t**3/6),
                     'K_exact_residual': abs(integrated_K-t**5/120)})

    for a in map(F, ('5/4', '3/2', '2', '3')):
        m, M, g, previous = F(1), F(2), F(1, 2), None
        for k in (1000, 2000, 4000, 8000, 16000):
            A, B, _ = matrix(a*a, m, M, g, F(k*k))
            _, _, D = eigenvalues(A, B, g)
            delta = float(A-B)
            weight = 2*float(g*g)/(D*(D+delta))
            reference = float(g*g/((a*a-1)**2*k**4))
            relative = abs(weight/reference-1)
            conditions = {'fast branch chi weight between zero and one': 0 < weight < 1,
                          'safe formula positive branch': delta > 0,
                          'asymptotic residue deviation': relative < ASYMPTOTIC_TOL}
            metrics = {'asymptotic_relative_deviation': relative}
            if previous is not None:
                ratio_deviation = abs((previous/weight)/16-1)
                conditions['doubling k^-4 ratio'] = ratio_deviation < ASYMPTOTIC_TOL
                metrics['doubling_ratio_relative_deviation'] = ratio_deviation
                diagnostics['residue_doubling_ratios'] += 1
            checks.case('fast_residue', f'{a}:{k}', locals_inputs(a=a, m=m, M=M, g=g, k=k),
                        conditions, metrics)
            previous = weight

    boundaries = [(2, 1, 1, 0), (2, 1, 1, 1), (2, 1, 1, -1),
                  (2, 0, 1, 0), (2, 1, 0, 0), (2, 0, 0, 0),
                  (1, 1, 1, F(1, 2)), (1, 1, 1, F(-1, 2))]
    boundary_reports = []
    for index, values in enumerate(boundaries):
        a, m, M, g = map(F, values)
        for k2 in (F(0), F(1)):
            conditions, metrics, details = mode_checks(
                a, m, M, g, [k2, F(0), F(0)], [F(2), F(-3)], [F(1), F(-2)],
                [F(1, 3), F(0), F(0)])
            if m == 0 or M == 0:
                conditions['zero mass requires zero mixing'] = g == 0
            if m > 0 and M > 0 and g*g == m*m*M*M:
                sign = F(1) if g > 0 else F(-1)
                null = (M, -sign*m)
                A0, B0, off0 = matrix(a*a, m, M, g, F(0))
                conditions['saturated mass null direction'] = (
                    A0*null[0]+off0*null[1] == 0 and off0*null[0]+B0*null[1] == 0)
                conditions['saturated square potential'] = (
                    hamiltonian(a*a, m, M, g, F(0), [F(2), F(-3)], [F(0), F(0)]) ==
                    (m*2+sign*M*(-3))**2/2)
            if g == 0:
                A, B, off = matrix(a*a, m, M, g, k2)
                hi, lo, _ = eigenvalues(A, B, off)
                conditions['g=0 decoupled spectrum'] = (
                    abs(hi-float(max(A, B))) <= FLOAT_TOL and abs(lo-float(min(A, B))) <= FLOAT_TOL)
                conditions['g=0 Born corrections vanish'] = leading_responses(g, F(3), F(5)) == (0, 0)
            if a == 1:
                A, B, _ = matrix(a*a, m, M, g, k2)
                conditions['coincident principal cone shifts mass matrix'] = A-k2 == m*m and B-k2 == M*M
            checks.case('mass_boundaries', f'{index}:{k2}', locals_inputs(a=a, m=m, M=M, g=g, k2=k2),
                        conditions, metrics)
            boundary_reports.append(encode({'a': a, 'm': m, 'M': M, 'g': g, 'k2': k2, **details}))

    amplitudes = list(map(F, (1, 2, 4, 8)))
    # Supercritical potential: minimize chi at fixed phi using the actual mass matrix.
    a2, m, M, g = F(4), F(1), F(2), F(3)
    A, B, off = matrix(a2, m, M, g, F(0))
    hi, lo, _ = eigenvalues(A, B, off)
    energies = [hamiltonian(a2, m, M, g, F(0), [n, -g*n/(M*M)], [F(0), F(0)]) for n in amplitudes]
    checks.case('pathological_controls', 'supercritical mixing', {'a2': a2, 'm': m, 'M': M, 'g': g, 'energies': energies},
                {'mass determinant negative': A*B-off*off < 0,
                 'negative homogeneous eigenvalue': lo < 0 < hi,
                 'negative scalable potential': all(e < 0 for e in energies),
                 'quadratic unbounded direction': all(energies[j] == energies[0]*amplitudes[j]**2 for j in range(4))})
    # Ghost energy is a kinetic-sign failure, independent of a negative potential.
    energies = [hamiltonian(F(4), F(1), F(1), F(0), F(0), [F(0), F(0)], [n, F(0)],
                            kinetic=(F(-1), F(1))) for n in amplitudes]
    healthy = [hamiltonian(F(4), F(1), F(1), F(0), F(0), [F(0), F(0)], [n, F(0)]) for n in amplitudes]
    checks.case('pathological_controls', 'negative kinetic ghost', {'energies': energies, 'healthy_energies': healthy},
                {'negative kinetic curvature': energies[0]*2 < 0,
                 'healthy kinetic counterpart positive': all(e > 0 for e in healthy),
                 'actual kinetic sign reversal': all(e == -h for e, h in zip(energies, healthy)),
                 'quadratic unbounded kinetic direction': all(energies[j] == energies[0]*amplitudes[j]**2 for j in range(4))})
    gradient_eigenvalues = [eigenvalues(*matrix(F(-1), F(1), F(1), F(0), F(k*k)))[1]
                           for k in (4, 8, 16)]
    gradient_energy = [hamiltonian(F(-1), F(1), F(1), F(0), F(k*k), [F(1), F(0)], [F(0), F(0)])
                       for k in (4, 8, 16)]
    checks.case('pathological_controls', 'negative gradient', {'k': [4, 8, 16], 'lambda_low': gradient_eigenvalues,
                                                            'energies': gradient_energy},
                {'negative high-k eigenvalues': all(v < 0 for v in gradient_eigenvalues),
                 'actual negative gradient energy': all(e < 0 for e in gradient_energy),
                 'growing high-k instability': gradient_eigenvalues[2] < gradient_eigenvalues[1] < gradient_eigenvalues[0]})
    A, B, off = matrix(F(4), F(0), F(0), F(1, 2), F(0))
    hi, lo, _ = eigenvalues(A, B, off)
    energy = hamiltonian(F(4), F(0), F(0), F(1, 2), F(0), [F(1), F(-1)], [F(0), F(0)])
    checks.case('pathological_controls', 'massless nonzero mixing', {'m': 0, 'M': 0, 'g': '1/2', 'energy': energy,
                                                                 'lambda_high': hi, 'lambda_low': lo},
                {'negative mass determinant': A*B-off*off < 0,
                 'negative homogeneous eigenvalue': lo < 0 < hi,
                 'negative potential direction': energy < 0})
    return diagnostics, boundary_reports


def locals_inputs(**kwargs):
    return kwargs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path, help='Explicit fresh JSON output path; overwrite refused.')
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        parser.error('Python 3.10 or later is required')
    if args.out.exists():
        parser.error('Output already exists; preserve it and choose a new --out path')
    checks, started = Checks(), timestamp()
    root = Path(__file__).resolve().parent
    registrations = {str(p.relative_to(root.parent)): sha(p) for p in (
        root.parent/'REGISTRATION.md', root.parent/'REGISTRATION_ORIGINAL.md',
        root/'IMPLEMENTATION_REGISTRATION.md', root/'IMPLEMENTATION_NOTES.md')}
    exception = None
    diagnostics, boundaries = {}, []
    try:
        diagnostics, boundaries = run(checks)
    except Exception as exc:
        exception = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    expected = {**REQUIRED, **EXTRA}
    count_checks = {name: checks.sections.get(name, {}).get('cases') == count
                    for name, count in expected.items()}
    count_checks['interior_grid_evaluations'] = diagnostics.get('interior_grid_evaluations') == 960
    count_checks['residue_doubling_ratios'] = diagnostics.get('residue_doubling_ratios') == 16
    all_passed = not exception and not checks.failures and all(count_checks.values())
    result = {'started_utc': started, 'completed_utc': timestamp(), 'all_passed': all_passed,
              'scope': 'Disclosed synthetic internal verification of an assumed coupled preferred-time scalar analogue.',
              'dependencies': {'python': platform.python_version(), 'libraries': 'Python standard library only'},
              'script_sha256': sha(Path(__file__).resolve()), 'registration_sha256': registrations,
              'parameters': {'mode_seed': 20261005, 'partial_fraction_seed': 20261006,
                             'spatial_integral_seed': 20261008, 'float_tolerance': FLOAT_TOL,
                             'asymptotic_relative_tolerance': ASYMPTOTIC_TOL,
                             'residual_scaling': 'Pre-run IMPLEMENTATION_REGISTRATION.md; raw and normalized maxima reported.',
                             'fast_residue_subset': {'a': ['5/4', '3/2', '2', '3'], 'm': 1, 'M': 2,
                                                     'g': '1/2', 'k': [1000, 2000, 4000, 8000, 16000]}},
              'expected_cases': expected, 'actual_count_checks': count_checks,
              'sections': checks.sections, 'diagnostics': diagnostics, 'boundary_cases': boundaries,
              'failed_cases': checks.failures, 'exception': exception,
              'qualifications': ['Massless nonzero mixing is unstable; H/K are Born coefficients only.',
                                  'Saturated homogeneous null modes can grow linearly; zero frequency is not strictly positive.',
                                  'Real-space kernels here are 1+1 dimensional distributions, not finite-energy pulses or 3+1 detector predictions.',
                                  'Continuum high-k modal suppression does not determine an EFT ultraviolet physical front.',
                                  'No apparatus cost, interaction with real matter, empirical FTL channel or external peer review is established.']}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    import json
    with args.out.open('x') as output:
        json.dump(encode(result), output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({'all_passed': all_passed, 'sections': checks.sections,
                      'failed_cases': len(checks.failures), 'exception': exception,
                      'script_sha256': result['script_sha256'], 'out': str(args.out)}, indent=2))
    return 0 if all_passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
