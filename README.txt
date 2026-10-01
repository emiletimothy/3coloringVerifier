Proposition 17: reproducible interval verification
================================================

Files
-----
verify.py               Standalone verifier (Python standard library only).
results.json            Output of the completed verification, including exact
                        rational margins, box counts, and a source SHA-256.
run.log                 Captured stdout from that execution.
verification_note.tex   Pasteable paragraph documenting the executed check.

Reproduce
---------
Use Python 3.10 or later, without the -O optimization flag:

    python3 verify.py --output my_results.json

Successful completion ends with:

    ALL REQUESTED CHECKS: PASS

An unresolved box, arithmetic exception, depth limit, or exhausted box budget
causes a nonzero exit and FAIL_OR_UNRESOLVED in the output. Such an outcome is
not evidence of a counterexample, but it is never reported as a certificate.
The entire recorded run took approximately 2.4 seconds on this machine.

Results of the recorded run
---------------------------
Initialization: PASS, 16,216 boxes visited, 8,140 final boxes.
  20 boxes were excluded by necessary feasibility constraints.
  148 boxes certified the first-lift threshold.
  7,972 boxes certified the second-lift threshold uniformly over successors.
  The minimum first-lift margin above T was > 2.16136e-5.
  The minimum second-lift margin above T was > 2.55026e-7.
  Thus every nonexcluded box certifies a threshold margin > 2e-7.

Drift: PASS, 646 boxes visited, 355 final boxes.
  15 boxes were excluded by necessary feasibility constraints.
  28 boxes passed the global Cauchy--Schwarz bound.
  312 boxes passed the convex supporting-line bound.
  The minimum certified lower bound on

      V * (M(q,L,z) - L - 1/1000)

  was > 2.83916e-5, uniformly over all admissible L and every
  c in [1/4, 250001/1000000]. In particular it exceeds 2e-5.

These are conservative lower bounds certified on subdivision boxes, not
estimates of the actual global minima. Exact fractions appear in results.json.
All original rectangles were covered. Both computations left zero unresolved
boxes. No random sampling is used to certify these statements.

Parameters and formulas
-----------------------
c_star = 250001/1000000, d = c/(1+c), D = 3(1+c), T = 53/25.
R_star = (4+3c)/D, S = R_star+1/10000.
The code uses the algebraically identical formula

    K(q) = min(5+3c, 4+3c + 3(1+c)*(1/10000 + rho(q))).

rho(q) is exactly the piecewise function in the manuscript:
    (1+2q)/(3-6q)                              if -1/2 <= q <= 0;
    (1-4q^2)/3                                if 0 < q < 1/16;
    min((1-4q^2)/3, (1-4q)/(2+4q))             if 1/16 <= q < 1/4;
    0                                         if q >= 1/4.
It is continuous and increases to q=0, then decreases. Evaluating endpoints
and zero (when present) therefore encloses it on each interval.

The lift is

    M(q,L,z) = ((1+2qz)L + 3sqrt(1-q^2)
                 - sqrt(K(q)-L^2)*sqrt(3d(1-q^2)-(2z+q)^2))
                / (2sqrt(1-q^2)*sqrt(1-z^2)).

Its domain has q,z in [-1/2,3/5], q+z <= 1/2,
(2z+q)^2 <= 3d(1-q^2), and 0 <= L <= sqrt(K(q)).

Arithmetic and coverage
-----------------------
All interval endpoints are integer multiples of 2^(-120). Addition,
subtraction, multiplication, division, squaring, and square roots use integer
arithmetic with outward rounding. Square-root endpoints use math.isqrt, an
integer square root. Every rational constant is enclosed outward too.

There are no floating-point acceptance or infeasibility tests. A small
floating-point bisection only proposes the rational ell used in the second-
lift dual bound. This proposal has no trusted mathematical role: ANY rational
ell is valid, and its resulting bound must pass an exact interval check.
Floating point is also used for human-readable progress and approximate
output margins. Exact fractions are authoritative.

Each rectangle starts as an 8-by-8 grid and unresolved boxes are bisected
along the larger relative width. All grid and bisection endpoints are exact
dyadic values. Children cover their parent, including the shared boundary.
The traversal must finish with an empty stack. Integer area accounting and
leaf-count identities are checked as additional consistency conditions.
Final-box records are hashed in deterministic traversal order; results.json
records that hash. The program regenerates the whole subdivision on each run.

The rho enclosure relies on the elementary monotonicity stated above.
The arithmetic self-test compares enclosure operations against Fraction
calculations at endpoints and interior points, including negative inputs.
This is a diagnostic in addition to the interval implementation, not a
substitute for reviewing the implementation and its mathematical reductions.

Initialization reduction
------------------------
Fix c=c_star and cover a in [0,d], z in [-1/2,3/5]. The actual starting
rectangle slightly enlarges the a endpoint to an outward dyadic bound on d.
Let b=(1+3a)/4 and

    eta(a) = (9-3a-6sqrt(c)*sqrt(d-a^2))/(4sqrt(1-b^2)).

Evaluate the first lift L1=M(b,eta(a),z) by interval arithmetic. Necessary
constraints can exclude an entire box only when their interval bounds prove
it has no feasible points. The verifier also supports immediate cap exits.

If the first threshold is not already certified above T, replace L1 by its
strictly positive dyadic lower bound l. The derivative of M with respect to
L is positive for L>=0 on continuing states, since 1+2qz >= 2/5.
It is therefore enough to bound every second lift starting at (z,l).

For general q,L, put s=sqrt(1-q^2), A=sqrt(3d)*s, W^2=K(q)-L^2.
Equation (13) in the PDF bounds H=V*(M(q,L,w)-T) below, for ANY rational ell,
by

    s^2 L + 3s + ell*q
      - A*sqrt((qL-ell)^2 + W^2)
      - 2sqrt(ell^2 + T^2 s^2).

This removes the second-successor variable w entirely. Since 0<V<=2,
a lower bound H>2*epsilon certifies M>T+epsilon. The code checks this with
epsilon=1/5000000. The bound is applied at q=z, L=l, exactly as in the paper.

When a radicand interval crosses zero, its lower endpoint is clamped to zero
ONLY under the corresponding feasibility/continuation constraint. This
intersects an enclosure with [0,infinity); it is valid for the feasible
subset. If the entire radicand enclosure is negative, that branch cannot
be used. In particular, this does not silently evaluate an infeasible lift.

Drift reduction
---------------
Cover q,z in [-1/2,3/5]. The entire c range is retained as one interval, so
there is no sampling or endpoint-only assumption in the treatment of c.
Put s=sqrt(1-q^2), B=sqrt(3d(1-q^2)-(2z+q)^2),
V=2s*sqrt(1-z^2), a0=1+2qz-V, and K=K(q). Then

    H(L)=V*(M(q,L,z)-L)=a0*L-B*sqrt(K-L^2)+3s.

For T <= L <= sqrt(K), the code uses the larger of two valid lower bounds:

    3s - sqrt(K)*sqrt(a0^2+B^2)

and, when K>T^2,

    H(T) + min(0,H'(T))*(sqrt(K)-T).

The second bound follows from convexity of H. The first has no division by
sqrt(K-T^2), so it also applies when K=T^2. No box is accepted merely because
a denominator is problematic. The code subtracts V/1000 and requires its
lower endpoint to exceed 1/50000 throughout the entire box.
