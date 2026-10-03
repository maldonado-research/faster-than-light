#!/usr/bin/env python3
"""Independent registered coupled-scalar checks, Python standard library only."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import platform
import random


BASE = Path(__file__).resolve().parent
ROOT_REGISTRATION = BASE.parent/'REGISTRATION.md'


def enc(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {str(k): enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def poly_value(coeff, x):
    return sum(c*x**i for i, c in enumerate(coeff))


def poly_derivative(coeff, order=1):
    for _ in range(order):
        coeff = [i*coeff[i] for i in range(1, len(coeff))]
    return coeff


def poly_integral(coeff, lo, hi):
    return sum(c*(hi**(i+1)-lo**(i+1))/Q(i+1) for i, c in enumerate(coeff))


def kernel_polys(a, t):
    d = a*a-1
    h_in = [t*t/(4*(a+1)), Q(0), -1/(4*a*(a+1))]
    h_out = [a*a*t*t/(4*a*d), -2*a*t/(4*a*d), 1/(4*a*d)]
    k_in = [a*(2*a+1)*t**4/(96*a*(a+1)**2), Q(0),
            -6*a*t*t/(96*a*(a+1)**2), Q(0), (a+2)/(96*a*(a+1)**2)]
    k_out = [Q(math.comb(4, i))*(a*t)**(4-i)*(-1)**i/(48*a*d*d) for i in range(5)]
    return h_in, h_out, k_in, k_out


def kernels(a, t, x):
    r = abs(x)
    if t <= 0 or r >= a*t:
        return Q(0), Q(0)
    hi, ho, ki, ko = kernel_polys(a, t)
    if r <= t:
        return poly_value(hi, r), poly_value(ki, r)
    return poly_value(ho, r), poly_value(ko, r)


def radial3_scaled(a, t, r):
    """Return 4*pi times formal 3D kernels, r>0, no division at origin."""
    if not r > 0:
        raise ValueError('radial dimension raising requires r>0')
    if t <= 0 or r >= a*t:
        return Q(0), Q(0)
    if r <= t:
        return 1/(a*(a+1)), (3*a*t*t-(a+2)*r*r)/(12*a*(a+1)**2)
    delta = a*t-r
    return delta/(a*(a*a-1)*r), delta**3/(6*a*(a*a-1)**2*r)


def dispersion_matrix(c_phi_squared, c_chi_squared, mu, nu, g, q):
    """Shared preferred positive-kinetic dispersion matrix, including sign controls."""
    return [[c_phi_squared*q+mu, g], [g, c_chi_squared*q+nu]]


def quadratic_energy(c_phi_squared, c_chi_squared, mu, nu, g,
                     phi, chi, phi_t, chi_t, grad_phi, grad_chi,
                     kinetic_phi=1, kinetic_chi=1):
    """Hamiltonian in velocity variables; preferred kinetic signs are parameters."""
    return (kinetic_phi*phi_t**2+kinetic_chi*chi_t**2
            +c_phi_squared*sum(x*x for x in grad_phi)
            +c_chi_squared*sum(x*x for x in grad_chi)
            +mu*phi**2+nu*chi**2+2*g*phi*chi)/2


def eigenvalues_from_matrix(K):
    A, B, g = K[0][0], K[1][1], K[0][1]
    D = math.hypot(A-B, 2*g)
    plus = (A+B+D)/2
    det = A*B-g*g
    minus = det/plus if plus else 0.0
    return A, B, D, plus, minus


def eigenvalues(a, mu, nu, g, q):
    return eigenvalues_from_matrix(dispersion_matrix(a*a, 1, mu, nu, g, q))


def stable_mass_criterion(mu, nu, g):
    return mu >= 0 and nu >= 0 and g*g <= mu*nu


def mixed_response_candidate(a, t, x, g):
    """Leading off-diagonal and same-chi mixed responses, not the unmixed term."""
    h, k = kernels(a, t, x)
    return -g*h, g*g*k


def mixed_radial_candidate(a, t, r, g):
    h3, k3 = radial3_scaled(a, t, r)
    return -g*h3, g*g*k3


def registered_candidate_checker(candidate):
    """The identical registered property gate evaluates healthy and injected candidates."""
    evaluated, failed = [], []

    def require(name, condition, parameters):
        entry = {'property': name, 'parameters': enc(parameters), 'passed': bool(condition)}
        evaluated.append(entry)
        if not condition:
            failed.append(entry)

    for mu, nu, g in [(Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(1,2)),
                      (Q(0),Q(1),Q(0)),(Q(1),Q(0),Q(0)),
                      (Q(1),Q(4),Q(2)),(Q(1),Q(4),Q(-2)),
                      (Q(1),Q(1),Q(2)),(Q(1),Q(1),Q(-2)),(Q(1),Q(1),Q(1,2))]:
        # Independent spectral behavior at k=0 checks the candidate mass classification.
        least = eigenvalues_from_matrix(dispersion_matrix(4.,1.,float(mu),float(nu),float(g),0.))[4]
        require('mass_classification_matches_low_k_spectrum',
                candidate['mass_stable'](mu,nu,g) == (least >= 0), (mu,nu,g))
    for a in (Q(3,2),Q(2),Q(3)):
        for t in (Q(1,4),Q(1)):
            g=Q(1,2)
            inter_r=(a+1)*t/2
            cross, same=candidate['mixed_response'](a,t,inter_r,g)
            neg_cross, neg_same=candidate['mixed_response'](a,t,inter_r,-g)
            twice_cross, twice_same=candidate['mixed_response'](a,t,inter_r,2*g)
            require('cross_is_odd_same_chi_is_even',
                    cross == -neg_cross and same == neg_same, (a,t,inter_r,g))
            require('same_chi_is_quadratic_in_g',
                    twice_cross == 2*cross and twice_same == 4*same, (a,t,inter_r,g))
            for x in ((a+1)*t,-(a+1)*t):
                require('response_zero_outside_fast_cone',
                        candidate['mixed_response'](a,t,x,g) == (0,0), (a,t,x,g))
            require('response_zero_before_source',
                    candidate['mixed_response'](a,-t,inter_r,g) == (0,0), (a,-t,inter_r,g))
            for fraction in (Q(1,4),Q(1,2),Q(3,4)):
                r=(1+(a-1)*fraction)*t
                hi,ho,ki,ko=kernel_polys(a,t)
                actual=candidate['mixed_radial'](a,t,r,g)
                expected=(-g*(-2*poly_value(poly_derivative(ho),r)/r),
                          g*g*(-2*poly_value(poly_derivative(ko),r)/r))
                require('radial_response_is_dimension_raised',actual == expected,(a,t,r,g))
            delta=(a-1)*t/4
            r1,r2=a*t-delta,a*t-2*delta
            y1=candidate['mixed_radial'](a,t,r1,g)[1]*r1
            y2=candidate['mixed_radial'](a,t,r2,g)[1]*r2
            require('radial_same_chi_front_power_is_cubic', y2 == 8*y1 and y1 > 0,
                    (a,t,r1,r2,g))
    return {'passed': not failed, 'constraint_count': len(evaluated),
            'failed_constraints': failed, 'evaluated_constraints': evaluated}


def matrix_multiply(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def retarded_series(K, t):
    total = [[t, 0.0], [0.0, t]]
    term = [[t, 0.0], [0.0, t]]
    converged = False
    for n in range(1, 81):
        product = matrix_multiply(K, term)
        factor = -t*t/((2*n)*(2*n+1))
        term = [[factor*product[i][j] for j in range(2)] for i in range(2)]
        total = [[total[i][j]+term[i][j] for j in range(2)] for i in range(2)]
        if max(abs(v) for row in term for v in row) < 1e-16:
            converged = True
            break
    return total, converged


def sinc_omega(lam, t):
    return t if lam == 0 else math.sin(t*math.sqrt(lam))/math.sqrt(lam)


def retarded_spectral(a, mu, nu, g, q, t):
    A, B, D, plus, minus = eigenvalues(a, mu, nu, g, q)
    if D == 0:
        f = sinc_omega(plus, t)
        return [[f, 0.0], [0.0, f]]
    fp, fm = sinc_omega(plus, t), sinc_omega(minus, t)
    return [[((A-minus)*fp+(plus-A)*fm)/D, g*(fp-fm)/D],
            [g*(fp-fm)/D, ((B-minus)*fp+(plus-B)*fm)/D]]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    if args.out.exists() or args.out.is_symlink():
        ap.error('--out already exists; refusing to overwrite an existing result or symlink')
    rng = random.Random(2026100305)
    counts, failures = {}, []

    def test(name, ok, detail=None):
        counts[name] = counts.get(name, 0)+1
        if not ok:
            failures.append({'check': name, 'detail': enc(detail)})

    max_residual = 0.0
    for cell in range(600):
        a = Q(rng.randint(5, 12), 4)
        m, M = Q(rng.randint(1, 8), 4), Q(rng.randint(1, 8), 4)
        mu, nu = m*m, M*M
        g = Q(rng.randint(-8, 8), 8)*m*M
        k = [Q(rng.randint(-8, 8), 4) for _ in range(3)]
        q = sum(v*v for v in k)
        phi, chi, pp, pc = [Q(rng.randint(-8, 8), 4) for _ in range(4)]
        fp = [Q(rng.randint(-8, 8), 4) for _ in range(3)]
        fc = [Q(rng.randint(-8, 8), 4) for _ in range(3)]
        V = (mu*phi*phi+nu*chi*chi+2*g*phi*chi)/2
        completion = ((m*phi+(g/m)*chi)**2+(nu-g*g/mu)*chi*chi)/2
        test('rational_mass_and_potential', mu*nu-g*g >= 0 and V == completion and V >= 0)
        e = quadratic_energy(a*a,1,mu,nu,g,phi,chi,pp,pc,fp,fc)
        v = Q(rng.randint(-9, 9), 10)
        charge = e+v*(pp*fp[0]+pc*fc[0])
        sq = ((pp+v*fp[0])**2+(pc+v*fc[0])**2
              +(a*a-v*v)*fp[0]**2+(1-v*v)*fc[0]**2
              +a*a*sum(z*z for z in fp[1:])+sum(z*z for z in fc[1:]))/2+V
        test('rational_charge_square', charge == sq and charge >= 0)
        flux = -a*a*pp*fp[0]-pc*fc[0]
        test('rational_finite_speed_flux', abs(flux) <= a*e)
        shifted_det = ((a*a-1)*q+mu)*nu-g*g
        test('rational_omega_at_least_light', shifted_det >= 0 and nu >= 0)
        A, B, D, lp, lm = eigenvalues(float(a), float(mu), float(nu), float(g), float(q))
        for lam in (lp, lm):
            residual = abs((A-lam)*(B-lam)-float(g*g))
            max_residual = max(max_residual, residual)
            test('floating_characteristic_residual', residual <= 1e-10,
                 {'cell': cell, 'residual': residual})
            test('floating_positive_observer_frequency', lam >= float(q)-1e-12
                 and (lam == 0 or math.sqrt(lam)-float(v*k[0]) > 0))
        test('rational_mode_trace_product',
             (a*a*q+mu)+(q+nu) == (a*a+1)*q+mu+nu
             and (a*a*q+mu)*(q+nu)-g*g
             == a*a*q*q+(mu+a*a*nu)*q+mu*nu-g*g)

    matrix_max_error = 0.0
    for cell in range(200):
        a = rng.randint(5, 12)/4
        m, M = rng.randint(1, 8)/4, rng.randint(1, 8)/4
        mu, nu = m*m, M*M
        g = rng.randint(-8, 8)/8*m*M
        q = rng.randint(0, 16)/4
        t = rng.randint(0, 16)/64
        direct, converged = retarded_series(dispersion_matrix(a*a,1,mu,nu,g,q), t)
        spectral = retarded_spectral(a, mu, nu, g, q, t)
        error = max(abs(direct[i][j]-spectral[i][j]) for i in range(2) for j in range(2))
        matrix_max_error = max(matrix_max_error, error)
        test('matrix_retarded_independent_routes', converged and error <= 1e-10,
             {'cell': cell, 'error': error})

    for cell in range(120):
        a = Q(rng.randint(5, 12), 4)
        t = Q(rng.randint(1, 12), 10)
        d = a*a-1
        hi, ho, ki, ko = kernel_polys(a, t)
        h_norm = 2*(poly_integral(hi, 0, t)+poly_integral(ho, t, a*t))
        k_norm = 2*(poly_integral(ki, 0, t)+poly_integral(ko, t, a*t))
        test('exact_1d_kernel_normalization', h_norm == t**3/6 and k_norm == t**5/120)
        test('exact_cone_matching', poly_value(hi, t) == poly_value(ho, t)
             and poly_value(ki, t) == poly_value(ko, t)
             and poly_value(poly_derivative(hi), t) == poly_value(poly_derivative(ho), t)
             and poly_value(poly_derivative(ki), t) == poly_value(poly_derivative(ko), t))
        test('exact_front_vanishing', all(poly_value(p, a*t) == 0
             and poly_value(poly_derivative(p), a*t) == 0 for p in (ho, ko)))
        # Integrate r²*(4*pi G3) exactly; outer factor 1/r cancels one r.
        h3in = [Q(0), Q(0), 1/(a*(a+1))]
        h3out = [Q(0), a*t/(a*d), -1/(a*d)]
        k3in = [Q(0), Q(0), 3*a*t*t/(12*a*(a+1)**2), Q(0),
                 -(a+2)/(12*a*(a+1)**2)]
        k3out = [Q(0)]+[Q(math.comb(3,i))*(a*t)**(3-i)*(-1)**i/(6*a*d*d) for i in range(4)]
        h3norm = poly_integral(h3in, 0, t)+poly_integral(h3out, t, a*t)
        k3norm = poly_integral(k3in, 0, t)+poly_integral(k3out, t, a*t)
        test('exact_3d_radial_normalization', h3norm == t**3/6 and k3norm == t**5/120)
        for ratio in (Q(0), Q(1,2), Q(1), (a+1)/2, a, a+1):
            r = ratio*t
            h, kk = kernels(a, t, r)
            test('exact_kernel_evenness_and_nonnegative', (h, kk) == kernels(a, t, -r)
                 and h >= 0 and kk >= 0)
            test('exact_negative_time_support', kernels(a, -t, r) == (0,0))
            if r >= a*t:
                test('exact_outside_fast_support', h == 0 and kk == 0)
            elif r > t:
                delta = a*t-r
                test('exact_intercone_coefficients', h == delta**2/(4*a*d)
                     and kk == delta**4/(48*a*d*d))
            if r < a*t:
                hp, kp = (hi, ki) if r <= t else (ho, ko)
                light2 = (t*t-r*r)/8 if r <= t else Q(0)
                test('exact_kernel_second_derivative',
                     poly_value(poly_derivative(kp, 2), r) == (h-light2)/d)
                if r > 0:
                    h3, k3 = radial3_scaled(a, t, r)
                    test('exact_dimension_raising', h3 == -2*poly_value(poly_derivative(hp),r)/r
                         and k3 == -2*poly_value(poly_derivative(kp),r)/r)
            if r > 0:
                test('exact_3d_even_support', radial3_scaled(a, -t, r) == (0,0))
            test('exact_g_parity', Q(1,2)**2*kk == Q(-1,2)**2*kk
                 and -Q(1,2)*h == -(-Q(-1,2)*h))

    # Independently derived equal-speed convolution coefficients.
    for i in range(30):
        t = Q(i+1, 10)
        r = Q(i % 6, 5)*t
        sigma = max(t*t-r*r, Q(0))
        h1, k1 = sigma/8, sigma*sigma/128
        # At a=1 D=(d_t²-d_x²); acting once on K returns H.
        # For polynomials in sigma: D f(sigma)=4(sigma f''+f').
        d_k = 4*(sigma*Q(2,128)+2*sigma/Q(128))
        test('coincident_cone_limit', d_k == h1 and h1 >= 0 and k1 >= 0)

    residue_max_error = 0.0
    for a in (1.5, 2.0, 3.0):
        for mu in (1.0, 4.0):
            for nu in (1.0, 4.0):
                for sign in (-1, 1):
                    g = sign*math.sqrt(mu*nu)/2
                    for k in (1000.0, 2000.0, 4000.0):
                        delta = (a*a-1)*k*k+mu-nu
                        D = math.hypot(delta, 2*g)
                        weight = 2*g*g/(D*(D+delta))
                        leading = g*g/((a*a-1)**2*k**4)
                        error = abs(weight/leading-1)
                        residue_max_error = max(residue_max_error, error)
                        test('fast_branch_residue_asymptotic', error <= 1e-4)

    certificates = []
    for a in (Q(3,2), Q(2), Q(3)):
        for g in (Q(-1,2), Q(1,2)):
            for t in (Q(1,100), Q(1,20), Q(1,10)):
                for fraction in (Q(1,4), Q(1,2), Q(3,4)):
                    r = (1+(a-1)*fraction)*t
                    _, k0 = kernels(a,t,r)
                    b = 1-t*t/4
                    lower = g*g*b**3*k0
                    # Uniform entire Volterra majorant; n=1..6 retained.
                    upper_prefix = sum(g**(2*n)*t**(4*n)/(2*a*math.factorial(4*n)) for n in range(1,7))
                    first_tail = g**14*t**28/(2*a*math.factorial(28))
                    ratio_bound = g*g*t**4/Q(29*30*31*32)
                    tail_upper = first_tail/(1-ratio_bound)
                    upper = upper_prefix+tail_upper
                    test('stable_full_volterra_positive_certificate', t <= 1 and 1-g*g >= 0
                         and t < r < a*t and lower > 0 and upper >= lower and 0 <= ratio_bound < 1)
                    certificates.append({'a':a,'m':1,'M':1,'g':g,'t':t,'r':r,
                                         'exact_full_response_lower_bound':lower,
                                         'uniform_majorant_upper_bound':upper,
                                         'majorant_tail_upper_bound':tail_upper,
                                         'bound_type':'analytic full stable 1+1 retarded response certificate'})

    for i in range(60):
        a = Q(rng.randint(5,12),4)
        mu = Q(rng.randint(1,8),2)
        nu = Q(rng.randint(1,8),2)
        g = Q(rng.randint(-3,3),10)
        p = Q(rng.randint(-10,10),100)*mu
        q = Q(rng.randint(0,8),4)
        w2 = a*a*q+p
        exact = -w2+q+nu-g*g/(mu-p)
        approximate = -w2+q+nu-g*g/mu-g*g*p/(mu*mu)
        remainder = -g*g*p*p/(mu*mu*(mu-p))
        Zt, Zx = 1+g*g/(mu*mu), 1+a*a*g*g/(mu*mu)
        test('exact_heavy_phi_remainder', exact-approximate == remainder and abs(p) <= mu/10)
        test('exact_heavy_phi_low_energy_coefficients',
             Zx/Zt >= 1 and ((Zx/Zt > 1) == (g != 0)))

    for a in (Q(3,2),Q(2),Q(3)):
        for j in range(1,17):
            # Rational Lorentz boosts and tensor transformations, not a
            # new propagating field state with the background held unchanged.
            p,qboost=Q(j),Q(5)
            v=(p*p-qboost*qboost)/(p*p+qboost*qboost)
            gamma=(p*p+qboost*qboost)/(2*p*qboost)
            A=gamma*gamma*(1-a*a*v*v)
            B=gamma*gamma*v*(a*a-1)
            C=gamma*gamma*(v*v-a*a)
            test('exact_boosted_principal_and_cauchy',A*C-B*B == -a*a
                 and gamma*gamma*(1-v*v) == 1
                 and (A > 0) == (abs(v)*a < 1)
                 and gamma*gamma*(1-v*v) == 1)  # chi kinetic coefficient
            if A < 0:
                # Leading transverse q_perp->infinity fast frequency square
                # is a²*q_perp²/A <0; algebraic mixing cannot change principal order.
                test('noncauchy_transverse_principal_growth',a*a/A < 0)

    boundaries = [(Q(0),Q(0),Q(0)), (Q(0),Q(1),Q(0)), (Q(1),Q(0),Q(0)),
                  (Q(1),Q(4),Q(2)), (Q(1),Q(4),Q(-2)), (Q(1),Q(4),Q(0))]
    for mu,nu,g in boundaries:
        A,B,D,lp,lm=eigenvalues(2.,float(mu),float(nu),float(g),0.)
        test('boundary_stability_and_zero_mode', lm >= 0 and ((lm == 0) == (mu*nu == g*g)))
        response,converged=retarded_series([[float(mu),float(g)],[float(g),float(nu)]],0.125)
        spectral=retarded_spectral(2.,float(mu),float(nu),float(g),0.,0.125)
        test('boundary_retarded_degeneracy', converged and max(abs(response[i][j]-spectral[i][j])
             for i in range(2) for j in range(2)) <= 1e-10)
    # Keep controls classified independently of healthy cases.
    control_outcomes = {}
    for name,mu,nu,g in [('massless_nonzero_g',0.,0.,0.5),
                         ('excessive_positive_g',1.,1.,2.),('excessive_negative_g',1.,1.,-2.)]:
        A,B,D,lp,lm=eigenvalues(2.,mu,nu,g,0.)
        rejected = lm < 0 and mu*nu-g*g < 0
        control_outcomes[name]={'type':'tachyonic_low_k','rejected':rejected}
        test('pathological_low_k_controls',rejected)
    bad_gradient_matrix=dispersion_matrix(-Q(4),Q(1),Q(1),Q(1),Q(1,2),Q(100))
    good_gradient_matrix=dispersion_matrix(Q(4),Q(1),Q(1),Q(1),Q(1,2),Q(100))
    bad_gradient_mode=eigenvalues_from_matrix(bad_gradient_matrix)[4]
    good_gradient_mode=eigenvalues_from_matrix(good_gradient_matrix)[4]
    gradient_state=(Q(0),Q(0),Q(0),Q(0),[Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(0)])
    bad_gradient_energy=quadratic_energy(-Q(4),Q(1),Q(1),Q(1),Q(1,2),*gradient_state)
    good_gradient_energy=quadratic_energy(Q(4),Q(1),Q(1),Q(1),Q(1,2),*gradient_state)
    control_outcomes['negative_gradient']={
        'type':'high_k_instability','rejected':bad_gradient_mode < 0 and good_gradient_mode > 0
        and bad_gradient_energy < 0 and good_gradient_energy > 0,
        'negative_gradient_matrix':bad_gradient_matrix,'negative_gradient_least_lambda':bad_gradient_mode,
        'healthy_least_lambda':good_gradient_mode,'negative_gradient_state_energy':bad_gradient_energy,
        'healthy_state_energy':good_gradient_energy}
    test('pathological_high_k_control',control_outcomes['negative_gradient']['rejected'])
    kinetic_state=(Q(0),Q(0),Q(1),Q(0),[Q(0)]*3,[Q(0)]*3)
    ghost_energy=quadratic_energy(Q(4),Q(1),Q(1),Q(1),Q(1,2),*kinetic_state,kinetic_phi=-1)
    healthy_energy=quadratic_energy(Q(4),Q(1),Q(1),Q(1),Q(1,2),*kinetic_state,kinetic_phi=1)
    control_outcomes['negative_preferred_kinetic']={
        'type':'ghost','rejected':ghost_energy < 0 and healthy_energy > 0,
        'negative_kinetic_state_energy':ghost_energy,'healthy_state_energy':healthy_energy}
    test('pathological_ghost_control',control_outcomes['negative_preferred_kinetic']['rejected'])
    healthy_candidate={'mass_stable':stable_mass_criterion,
                       'mixed_response':mixed_response_candidate,'mixed_radial':mixed_radial_candidate}
    candidate_gate=registered_candidate_checker(healthy_candidate)
    test('registered_candidate_healthy',candidate_gate['passed'],candidate_gate['failed_constraints'])

    def wrong_mass_criterion(mu,nu,g):
        return stable_mass_criterion(mu,nu,g) or (mu == 0 and nu == 0)

    def wrong_odd_same_chi(a,t,x,g):
        h,k=kernels(a,t,x)
        return -g*h,g*k

    def wrong_outside_response(a,t,x,g):
        if t > 0 and abs(x) > a*t:
            return -g,g*g
        return mixed_response_candidate(a,t,x,g)

    def wrong_copied_1d_radial(a,t,r,g):
        h3,_=radial3_scaled(a,t,r)
        _,k1=kernels(a,t,r)
        return -g*h3,g*g*k1

    substitutions={
        'massless_mixing_called_stable':('mass_stable',wrong_mass_criterion),
        'same_chi_wrong_odd_g':('mixed_response',wrong_odd_same_chi),
        'outside_fast_nonzero':('mixed_response',wrong_outside_response),
        'copied_1d_front_power_to_3d':('mixed_radial',wrong_copied_1d_radial)}
    mutation_controls={}
    for name,(entry,replacement) in substitutions.items():
        injected=dict(healthy_candidate)
        injected[entry]=replacement
        gate=registered_candidate_checker(injected)
        mutation_controls[name]={'injected_entry':entry,'injected_function':replacement.__name__,
                                 'rejected':not gate['passed'],'gate':gate}
        test('registered_mutant_rejection',not gate['passed'],name)

    result={'status':'PASS' if not failures else 'FAIL',
            'evidence_tier':'internal registered synthetic mathematical verification',
            'seed':2026100305,'python':platform.python_version(),'dependencies':'Python standard library only',
            'script_sha256':sha(Path(__file__)), 'independent_registration_sha256':sha(BASE/'REGISTRATION.md'),
            'root_registration_sha256_at_run':sha(ROOT_REGISTRATION),
            'root_implementation_read_before_run':False,
            'counts':counts,'total_checks':sum(counts.values()),'failures':failures,
            'maximum_characteristic_absolute_residual':max_residual,
            'maximum_retarded_matrix_absolute_error':matrix_max_error,
            'maximum_fast_residue_relative_error':residue_max_error,
            'stable_intercone_full_response_certificates':certificates,
            'pathological_controls':control_outcomes,'healthy_candidate_gate':candidate_gate,
            'mutants_rejected':mutation_controls,
            'correction_history':{
                'initial_results_sha256':sha(BASE/'INITIAL_RESULTS.json'),
                'initial_script_sha256':sha(BASE/'INITIAL_CHECK.py'),
                'initial_count':10414,
                'initial_weak_mutant_comparisons':'Four truth-versus-wrong-value comparisons did not inject candidates and are not certified mutant rejection.',
                'initial_gradient_kinetic_controls':'Literal negative expectations were replaced by shared parameterized dispersion/energy helper evaluations.',
                'final_change':'Same four named defects injected into candidate functions and fed through identical 63-constraint gate; healthy gate added.',
                'count_delta':sum(counts.values())-10414,
                'not_blind':'Post-run amendment after skeptical review; initial artifacts preserved.'},
            'scope':['1+1 massive stable exact-response lower bounds come from a convergent positive Volterra series.',
                     '3D radial formulas require r>0 and are mathematical kernels, not experiment predictions.',
                     'Massless nonzero mixing is unstable; its kernels serve as perturbation coefficients only.',
                     'Full-continuum front statements do not settle unknown EFT ultraviolet completion.',
                     'No matter/photon detector, observation, physical FTL, novelty or external review claim.']}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as output:
        output.write(json.dumps(enc(result),indent=2)+'\n')
    print(f"{result['status']}: {result['total_checks']} checks, {len(certificates)} full-response certificates; out={args.out}")
    return 1 if failures else 0


if __name__=='__main__':
    raise SystemExit(main())
