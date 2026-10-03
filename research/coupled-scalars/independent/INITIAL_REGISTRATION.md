# Independent coupled-scalar verification registration

2026-10-03 UTC, prepared after the root registration gate and before this
reviewer's decisive calculations. The root registration was reported bound
at 00:37:27 UTC with SHA-256
`b86cb57adec77f12a7a76d398fdafcdf1e65a5ae7f81776148b69cee7f190718`.
Later disclosed root amendments are separately hashed at execution.

No root implementation has been read. Before this registration, this reviewer
derived the mass-matrix criterion, spectra, preferred/boosted energy bounds,
massless piecewise convolution polynomials, 3D radial-derivative identity,
front powers, and small-time massive positivity certificate. The root supplied
the suggested low-energy heavy-phi elimination coefficients, which were then
independently algebraically verified. This is disclosed internal verification,
not blind preregistration or external review. No observational input is used.

The model and scope follow the parent registration: a>1, m,M>=0, real g,
ordinary c=1, constant global preferred time, positive preferred kinetic and
gradient terms, local quadratic mixing, mathematical source J_chi and measured
field chi. Neither scalar is identified with actual detector matter or photons.
Retarded kernels are distributions; physical production/detection, source work,
apparatus, quantization/gravity and UV completion remain unsupplied. The 1D
and formal r>0 3D mathematical kernels are explicitly distinguished.

Planned verification uses Python >=3.10, standard library only, one process,
seed 2026100305. All outputs stay in this independent directory or explicit
temporary --out locations. The implementation must require --out, hash-bind
its own code and both registrations, save failures, and return nonzero on any
failed criterion. Exact Fraction checks require equality. Float matrix/spectral
comparisons use absolute tolerance 1e-10; no relative claim is made at a zero
response. No parent script or output supplies a test oracle.

1. Generate 600 deterministic rational stable cases: rational a in [5/4,3],
   positive rational m,M in [1/4,2], g=p*m*M with p in [-1,1], rational spatial
   momentum components and arbitrary rational fields/momenta. Check exact
   mass determinant, completed-square energy, stable eigenvalue lower bound
   via K-k²I PSD, ordinary-observer charge square completion, source power/flux
   bounds and g parity. Floating eigenvalues use a cancellation-safe determinant
   quotient and must obey the characteristic equation to 1e-10.
2. Generate 200 bounded stable Fourier/time cases, t<=1/4 and k²<=4. Compare
   the retarded matrix from spectral projectors to an independently evaluated
   matrix power series sin(t sqrt(K))/sqrt(K), using 80 terms maximum and
   convergence threshold 1e-16. Treat zero/degenerate eigenvalues explicitly.
3. Use 120 exact rational kernel parameter cells and six radial positions
   r/t={0,1/2,1,(1+a)/2,a,a+1}. Independently implemented piecewise H=Ga*G1
   and K=G1*Ga*G1 must vanish outside the fast cone, match the intercone
   quadratic/quartic forms, obey the spatial second-derivative identity
   K''=(H-G1*G1)/(a²-1), and normalize to t³/6 and t⁵/120 by exact polynomial
   integration across both cone regions. Verify continuity at r=t and vanishing
   value/slope at r=a*t. Both signs of x and g are checked.
4. In the same cells, verify the rotationally invariant dimension-raising
   relation G3=-(2*pi*r)^(-1)*partial_r G1 for r>0 using scaled 4*pi G3
   polynomial/rational formulas. Check 3D spatial normalization by exact
   radial integration and cubic/linear intercone powers for chi-chi/cross
   coefficients. These are mathematical kernel checks, not 3D experimental
   predictions. r=0 is excluded from the radial formula; distributional
   boundaries are derived rather than divided by zero. Test 30 coincident-
   cone a=1 points with independently derived sigma/8 and sigma²/128 limits.
5. Evaluate 72 large-k fast-branch chi spectral weights: a in {3/2,2,3},
   m²,M² in {1,4}, g=+/-m*M/2, k in {1000,2000,4000}. Use cancellation-safe
   residue 2g²/[D(D+Delta)] and require agreement with
   g²/[(a²-1)²k⁴] to relative 1e-4. This checks modal suppression, not true
   physical front or detector sensitivity.
6. Certify 54 stable nonzero-mixing intercone cases exactly: a in {3/2,2,3},
   m=M=1, g=+/-1/2, t in {1/100,1/20,1/10}, and
   r/t=1+(a-1)*{1/4,1/2,3/4}. For free 1D massive kernels and m*t<=1,
   (1-m²t²/4)G0<=Gm<=G0. Retarded convolution preserves positivity; all
   chi-to-chi Volterra terms carry g^(2n). Consequently the **exact full stable
   response**, not just an unstable massless Born coefficient, obeys
   G_chichi >= g²(1-t²/4)^3 K_massless >0 in the intercone. Convergence is
   certified by the majorant sum g^(2n)t^(4n)/[2a(4n)!], n>=1. Save each
   rational positive lower bound and an explicit finite majorant/tail bound.
   This is an analytic certificate whose parameter identities are computed,
   not a numerical measurement of the full response.
7. Check 60 exact rational heavy-phi elimination symbols with m²>0 and
   |p|<=m²/10, p=omega²-a²k². The difference between the exact eliminated
   symbol and its first derivative approximation must equal
   -g²p²/[m⁴(m²-p)]. Check Zt=1+g²/m⁴, Zx=1+a²g²/m⁴ and their speed ratio;
   label the truncation low-energy and do not infer a global front or ordinary
   rods/clocks from it.
8. Keep separately classified boundaries/controls: g=0, m=0, M=0, both
   masses zero, positive/negative saturation, a=1, k=0 zero eigenvalue,
   massless g!=0 tachyon, excessive mixing tachyon, negative gradient high-k
   instability and reversed preferred kinetic sign ghost. Mutants claiming
   massless nonzero mixing is stable, using an odd-in-g same-chi response,
   retaining support outside fast cone, or copying a 1D quartic onset into 3D
   must be rejected. Boundary cells are explicit, not silently skipped.

The analytic audit must explain that stable eigen-covectors have omega>=|k|,
ordinary observer translation charges remain nonnegative, common boosted
Cauchy slices require |v|a<1, and negative fast kinetic coefficient on a
non-Cauchy slice alone is not a preferred-frame ghost. Mixing changes chi
dispersion and makes its retarded support reach the fast cone; calling chi
unchanged ordinary matter would be misleading. Two principal characteristic
cones share one outer domain-of-dependence cone and one preferred temporal
orientation. Separate momentum-dependent branch projectors are nonlocal;
exterior cancellation concerns the complete retarded response.

Pre-run addendum, before the first decisive execution: include 48 exact
rational boosted principal-matrix cells, a in {3/2,2,3}, Pythagorean boosts
from p=1..16,q=5. Verify the fast 2x2 determinant, ordinary chi kinetic
coefficient 1, and common-Cauchy sign criterion |v|a<1. For each A_phi<0
cell verify negative leading transverse frequency square a²/A_phi; no
universal 1+1 wrong-slicing ill-posedness claim is made. Code is still unrun.
