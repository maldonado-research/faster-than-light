# Independent 1+1D retarded-kernel verification

Recorded 2026-10-03 UTC before this implementation's first decisive run. Read the root registration at ../REGISTRATION.md first. The kernel formulas, triangle transformation, Bessel bounds, Volterra convergence argument, and front coefficients were derived analytically before this record. This is disclosed internal verification, not blind discovery preregistration, external replication, physical detector modeling, or a novelty claim.

Scope is the explicitly assumed fixed-preferred-time two-scalar system. Source and response are the chi field only. All retarded impulse kernels and quadrature here are 1+1 dimensional. The source is idealized prescribed forcing; no finite apparatus energy, production mechanism, actual ordinary-matter identity, calibrated detector observable, UV completion, quantum/gravitational completion, or observation is supplied. No Git or Zenodo changes are part of this run.

Write scripts and results only in this kernel directory. Python >=3.10 and standard library only. No Fourier cutoff or numerical transform is used to infer support.

## Selected cases and fixed criteria

1. Evaluate 36 exact rational triangle moments, p,q=0,...,5. Compare polynomial integration after v=(1-u)z against p!q!/(p+q+2)!, including area 1/2 and integral(uv)=1/24. Exact equality is required.
2. Evaluate 24 independent rational transform cases: a in {5/4,3/2,2,3}, Laplace parameter s in {1/2,1,2}, Fourier k in {0,1}. Compare each original H and K denominator with its partial-fraction decomposition: 48 exact identities. Exact equality is required.
3. Use tensor composite Simpson resolutions 64 and 128 for nine stable massive cells: a=2,m=M=1,g=1/2, t in {1/100,1/20,1/10}, r/t in {5/4,3/2,7/4}. Calculate both H and K through the bounded intercone triangle integral, independently of root's implementation. Compare the normalized integral values at the two resolutions with relative tolerance 1e-9. Require their normalized ratios within the analytic [1-E,1] bounds, allowing numerical slack 2e-12.
4. Add five preselected near-front cells with a=2,m=M=1,g=1/2,t=1/10 and Delta/t in {1/2,1/4,1/8,1/16,1/32}. Use the same two resolutions and tolerances. Normalized H/H0 and K/K0 must increase toward 1. Four successive half-gap comparisons must recover quadratic H and quartic K scaling within relative 1e-2. At the final cell, (1-H/H0)/Delta and (1-K/K0)/Delta must agree with m^2*t/(6a) and m^2*t/(10a), respectively, within relative 1e-2.
5. Thus the massive/front calculation comprises 14 configured cell evaluations, 28 quadrature grids and 56 scalar integral values. There are 13 distinct parameter cells: the first front point repeats one of the prescribed stable cells. Separately evaluate two massless grids at resolutions 32 and 64, requiring the H and K triangle moments to be 1/2 and 1/24 within absolute 1e-13: four normalization checks.
6. Use twelve-term entire power series for J0(sqrt(z2)) and 2J1(sqrt(z2))/sqrt(z2), with the latter defined as 1 at z2=0. Require 0<=z2<=1. Alternating decreasing terms make the first omitted term a rigorous truncation bound, at most 1/[4^12*(12!)^2] for J0 and 1/[4^12*12!*13!] for the second factor. These analytic series bounds are distinct from floating roundoff and Simpson resolution agreement.
7. For all 14 massive/front cells, evaluate exact rational small-time chi lower bounds and convergent higher-order upper bounds. Require g^2<=m^2*M^2, m*t<=1, M*t<=1, a valid intercone point, strictly positive lower bound, and the higher-order geometric ratio strictly below 1. Quadrature is not the support proof: positivity and convergence of the even-mixing Volterra series supply that proof.
8. Classify three controls exactly: a=2,m=M=1,g=1/2 is stable; m=M=0,g=+/-1/2 has a negative k=0 eigenvalue and is unstable. The latter supplies only formal massless perturbative coefficients, never a stable full-model example.

Expected counts: 36 exact moments; 24 transform cases/48 identities; 14 massive/front cells/28 grids/56 scalar values; 2 massless grids/4 scalar normalization checks; 14 analytic full-response-bound cells; 8 adjacent normalized-ratio monotonicity checks; 8 power-scaling checks; 2 final slope checks; 3 stability controls.

Any failed assertion is recorded as failure and produces a nonzero exit status. Zero checks, exceptions, or skipped decisive calculations do not establish success. Preserve this file and append a dated amendment if criteria change after results are inspected. Record source/registration hashes, dependencies, parameters, counts, resolution changes, bounds, failures, and current output timestamp.

## Analytic premises

Outside the slow cone, Delta=a*t-r and delta=a^2-1. With u,v>=0,u+v<=1, set S=4*Delta^2*u*v/delta and

F=(Delta/a)*(1-u-v)*[2*t-(Delta/a)*(1+((a+1)/(a-1))*u+((a-1)/(a+1))*v)].

H=Delta^2/(2*a*delta) integral(J0(m*sqrt(F))*J0(M*sqrt(S)) du dv).

K=Delta^4/(2*a*delta^2) integral(u*v*J0(m*sqrt(F))*[2*J1(M*sqrt(S))/(M*sqrt(S))] du dv).

Their massless values are H0=Delta^2/(4*a*delta) and K0=Delta^4/(48*a*delta^2).

Relative errors from the massless coefficients are bounded by E_H=m^2*t*Delta/(2*a)+M^2*Delta^2/(4*delta), and E_K=m^2*t*Delta/(2*a)+M^2*Delta^2/(8*delta).

The locally uniformly convergent full chi series is sum(n>=0) g^(2n) G_1,M^{*(n+1)}*G_a,m^{*n}. Each coefficient is pointwise bounded in absolute value by t^(4n)/[2*(4n)!]. At m*t,M*t<=1 all factors and all even terms are nonnegative. In the intercone the n=0 term vanishes, giving

G_chichi >= g^2*(1-m^2*t^2/4)*(1-M^2*t^2/4)^2*K0 > 0.

A sharper intercone absolute coefficient bound is

K_n <= t^(n-1)*Delta^(3n+1)/[2^n*a^n*delta^(n+1)*(n-1)!*(3n+1)!], n>=1.

The n=1 value is exactly K0. For the higher-order remainder, let B2=g^4*t*Delta^7/[20160*a^2*delta^3] and eta=g^2*t*Delta^3/[2880*a*delta]. Since successive bounds decrease at least geometrically with ratio eta from n=2 onward, the remainder is <=B2/(1-eta) when eta<1. This supports a full-response upper bound g^2*K0+B2/(1-eta).

At fixed t,a,m,M, the leading normalized front slopes are 1-H/H0=m^2*t*Delta/(6a)+O(Delta^2) and 1-K/K0=m^2*t*Delta/(10a)+O(Delta^2). The full intercone chi response has leading g^2*Delta^4/[48*a*delta^2]; higher mixing is smaller by O(g^2*t*Delta^3). None of these bounds establishes observability.

## Source-only output portability amendment

Recorded 2026-10-03T00:57:07.456123+00:00 after inspection of the first successful run and before the final rerun. The original directory-only output guard prevented the public reproduction command from writing to a temporary directory. Root authorized this source-only correction: an explicit --out may now name any writable fresh file; existing files, symlinks, canonical results and source/registration records are refused. Creation uses exclusive mode. No formula, parameter, check count, quadrature resolution, mathematical bound, or acceptance tolerance changes. The earlier restriction on writing within this kernel directory applies to the development task rather than the portable script's reproduction output.

The initial script and output are retained byte-for-byte as INITIAL_REPRODUCE.py and INITIAL_RESULTS.json. Their SHA-256 values are, respectively:

- `1830892f8b9e3c0c6b1a449ba5224cfe34459e4e60e049931bfa148dbe194388`
- `d4dd2420ef4b47260281858f158b8e4aea4f09a1e6585439b89fbc151fecfb3c`

The registration before this amendment had SHA-256 `59d2de44809f56a0e8530cc7cacfe30877b75a867fc7d2a3dc752f5cff0fbee5`. The amended script at this point has SHA-256 `4bdd99418bdd372703f03ea275b7f2d840ebeeb61469ae6e1403c323f1d3d5f4`. One final decisive run will write to a fresh explicit /tmp path, bind the amended registration and source hashes, and be compared with the initial mathematical payload. Only after inspection will its bytes be explicitly promoted to canonical RESULTS.json. The preserved initial records are historical development evidence; their hashes do not identify the amended script.

Pre-rerun metadata correction, disclosed after review: the 14 configured massive/front evaluations contain **12 distinct parameter cells**, as the implemented result already reports. At t=1/10, the front points Delta/t=1/2 and 1/4 repeat the prescribed r/t=3/2 and 7/4 stable cells. The original registration's statement of 13 distinct cells and only one repeated point is retained above as historical text and superseded here. The actual 14 evaluations, 28 grids, 56 scalar integrals, checks, parameters and tolerances are unchanged. This corrects a descriptive registration count; no decisive case is added, removed or rerun selectively.
