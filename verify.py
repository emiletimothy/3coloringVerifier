#!/usr/bin/env python3
"""Exact dyadic interval verification of Proposition 17 ("finite inequalities").

Python 3.10+ standard library only. Run: python3 verify.py --output results.json
All acceptance/rejection decisions use integers at scale 2**120. Floating point
is used ONLY to suggest a rational dual witness; every such witness is checked
by exact interval arithmetic. An unresolved box is an error, never a pass.
See README.txt for the mathematical reductions and scope of the certificate.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import time

BITS = 120
Q = 1 << BITS


def ceildiv(a, b):
    assert b > 0
    return -((-a) // b)


class I:
    """Closed interval [lo/Q, hi/Q]; all endpoints are integers."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi
        assert isinstance(self.lo, int) and isinstance(self.hi, int)
        assert self.lo <= self.hi

    @classmethod
    def rat(cls, n, d=1):
        assert d > 0
        return cls((n*Q)//d, ceildiv(n*Q, d))

    @staticmethod
    def coerce(x):
        return I.rat(x) if isinstance(x, int) else x

    def __add__(self, other):
        other = self.coerce(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, int):
            return I(self.lo*other, self.hi*other) if other >= 0 else I(self.hi*other, self.lo*other)
        vals = (self.lo*other.lo, self.lo*other.hi,
                self.hi*other.lo, self.hi*other.hi)
        return I(min(vals)//Q, ceildiv(max(vals), Q))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval division through zero")
        if other.hi < 0:
            return (-self)/(-other)
        pairs = [(a*Q, b) for a in (self.lo, self.hi)
                 for b in (other.lo, other.hi)]
        return I(min(a//b for a, b in pairs),
                 max(ceildiv(a, b) for a, b in pairs))

    def sq(self):
        upper = max(self.lo*self.lo, self.hi*self.hi)
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        return I(lower//Q, ceildiv(upper, Q))

    def sqrt(self):
        assert self.lo >= 0, "negative square-root interval"
        lo = math.isqrt(self.lo*Q)
        hi = math.isqrt(self.hi*Q)
        if hi*hi < self.hi*Q:
            hi += 1
        return I(lo, hi)

    def nonnegative(self):
        """Intersect with [0,infinity); used ONLY under feasibility constraints."""
        assert self.hi >= 0
        return I(max(0, self.lo), self.hi)

    def mid_float(self):
        return (self.lo+self.hi)/(2*Q)


def imin(x, y):
    x, y = I.coerce(x), I.coerce(y)
    return I(min(x.lo, y.lo), min(x.hi, y.hi))


def rho_point(x):
    """Enclose rho at the exact dyadic point x/Q."""
    q = I(x)
    if x <= 0:
        return (1+2*q)/(3-6*q)
    if 16*x < Q:
        return (1-4*q.sq())/3
    if 4*x < Q:
        return imin((1-4*q.sq())/3, (1-4*q)/(2+4*q))
    return I(0)


def rho(q):
    # rho is continuous, increasing on [-1/2,0], decreasing on [0,3/5].
    # Both positive-branch candidates decrease; at 1/16 the minimum agrees
    # with the earlier branch. Include the only maximum, at zero.
    ends = [rho_point(q.lo), rho_point(q.hi)]
    if q.lo <= 0 <= q.hi:
        ends.append(I.rat(1, 3))
    return I(min(x.lo for x in ends), max(x.hi for x in ends))


CSTAR = I.rat(250001, 1000000)
CRANGE = I(I.rat(1, 4).lo, CSTAR.hi)
T = I.rat(53, 25)
INIT_EPS = I.rat(1, 5000000)  # 2e-7 threshold margin
DRIFT_EPS = I.rat(1, 50000)   # 2e-5 numerator margin beyond V/1000
DELTA = I.rat(1, 1000)


def cap(q, c):
    # D*S = 4+3c + D/10000, with D=3(1+c).
    return imin(5+3*c, 4+3*c+3*(1+c)*(I.rat(1, 10000)+rho(q)))


def eta(a, c, d):
    b = (1+3*a)/4
    return (9-3*a-6*c.sqrt()*(d-a.sq()).sqrt())/(4*(1-b.sq()).sqrt())


def suggest_ell(q, L, c):
    """Untrusted floating-point search for a rational dual witness."""
    qf, lf, cf = q.mid_float(), L.mid_float(), c.mid_float()
    s2 = 1-qf*qf
    A = math.sqrt(3*cf/(1+cf)*s2)
    W2 = max(0., cap(I((q.lo+q.hi)//2), c).mid_float()-lf*lf)
    lo, hi = -16., 16.
    for _ in range(40):
        ell = (lo+hi)/2
        diff = qf*lf-ell
        r = math.sqrt(diff*diff+W2)
        deriv = qf + (A*diff/r if r else 0.) - 2*ell/math.sqrt(ell*ell+(53/25)**2*s2)
        if deriv > 0:
            lo = ell
        else:
            hi = ell
    return I.rat(round(((lo+hi)/2)*(1 << 20)), 1 << 20)


def second_dual(q, L, c, d, ell):
    s2 = 1-q.sq()
    s = s2.sqrt()
    w2 = (cap(q, c)-L.sq()).nonnegative()
    A = (3*d*s2).sqrt()
    return s2*L+3*s+ell*q-A*((q*L-ell).sq()+w2).sqrt()-2*(ell.sq()+T.sq()*s2).sqrt()


def initialization(a, z):
    c = CSTAR
    d = c/(1+c)
    b = (1+3*a)/4
    # Discard only boxes certified disjoint from a necessary constraint.
    if (b+z).lo > I.rat(1, 2).hi:
        return "infeasible", None
    B2 = 3*d*(1-b.sq())-(2*z+b).sq()
    if B2.hi < 0:
        return "infeasible", None
    L0 = eta(a, c, d)
    K0 = cap(b, c)
    if (L0.sq()-K0).lo > INIT_EPS.hi:
        return "seed_cap", (L0.sq()-K0).lo
    if (L0-T).lo > INIT_EPS.hi:
        return "seed_threshold", (L0-T).lo
    W2 = K0-L0.sq()
    if W2.hi < 0:
        # Too small a terminal margin for our requested quantitative bound.
        return None, None
    s = (1-b.sq()).sqrt()
    V = 2*s*(1-z.sq()).sqrt()
    L1 = ((1+2*b*z)*L0+3*s-W2.nonnegative().sqrt()*B2.nonnegative().sqrt())/V
    if (L1-T).lo > INIT_EPS.hi:
        return "first_threshold", (L1-T).lo
    if L1.lo <= 0:
        return None, None
    # M increases in L>=0 on continuing states; replace L1 by exact lower end.
    low = I(L1.lo)
    if (low.sq()-cap(z, c)).lo > INIT_EPS.hi:
        return "first_cap", (low.sq()-cap(z, c)).lo
    if (cap(z, c)-low.sq()).hi < 0:
        return None, None
    ell = suggest_ell(z, low, c)
    bound = second_dual(z, low, c, d, ell)
    # H=V*(M-T), and 0<V<=2. Therefore H>2*eps gives M>T+eps.
    if bound.lo > (2*INIT_EPS).hi:
        return "second_threshold", bound.lo//2
    return None, None


def drift(q, z):
    c = CRANGE  # A single interval covers every c in [1/4,1/4+1e-6].
    d = c/(1+c)
    if (q+z).lo > I.rat(1, 2).hi:
        return "infeasible", None
    s2 = 1-q.sq()
    B2 = 3*d*s2-(2*z+q).sq()
    if B2.hi < 0:
        return "infeasible", None
    K = cap(q, c)
    W2 = K-T.sq()
    if W2.hi < 0:
        return "no_continuing_state", None
    s = s2.sqrt()
    B = B2.nonnegative().sqrt()
    V = 2*s*(1-z.sq()).sqrt()
    a0 = 1+2*q*z-V
    global_bound = 3*s-K.sqrt()*(a0.sq()+B2.nonnegative()).sqrt()-V*DELTA
    best, kind = global_bound.lo, "global"
    if W2.lo > 0:
        W = W2.sqrt()
        HT = a0*T-B*W+3*s
        HpT = a0+B*T/W
        tangent = HT+imin(0, HpT)*(K.sqrt()-T)-V*DELTA
        if tangent.lo > best:
            best, kind = tangent.lo, "tangent"
    # The global bound has no singular denominator and remains valid at K=T^2.
    if best > DRIFT_EPS.hi:
        return kind, best
    return None, None


def exact_str(scaled):
    return str(Fraction(scaled, Q))


def grid(lo, hi, n=8):
    pts = [lo+(hi-lo)*k//n for k in range(n+1)]
    return [I(pts[k], pts[k+1]) for k in range(n)]


def run_cover(name, check, xr, yr, max_boxes, max_depth):
    # Closed child intervals share a boundary; their union equals the parent.
    roots = [(x, y, 0) for x in grid(*xr) for y in grid(*yr)]
    stack = list(reversed(roots))
    counts = Counter()
    margins = {}
    digest = hashlib.sha256()
    visited = split = depth_max = 0
    area = 0
    initial_area = (xr[1]-xr[0])*(yr[1]-yr[0])
    start = last = time.monotonic()
    while stack:
        x, y, depth = stack.pop()
        visited += 1
        depth_max = max(depth_max, depth)
        reason, margin = check(x, y)
        if reason is not None:
            counts[reason] += 1
            area += (x.hi-x.lo)*(y.hi-y.lo)
            if margin is not None:
                margins[reason] = min(margin, margins.get(reason, margin))
            record = [x.lo, x.hi, y.lo, y.hi, reason, margin]
            digest.update((json.dumps(record, separators=(",", ":"))+"\n").encode())
        else:
            if visited >= max_boxes or depth >= max_depth:
                raise RuntimeError(f"UNRESOLVED {name} box at depth {depth}: "
                                   f"{exact_str(x.lo)}..{exact_str(x.hi)}, "
                                   f"{exact_str(y.lo)}..{exact_str(y.hi)}")
            split += 1
            # Bisect the larger relative width using integer comparisons.
            if (x.hi-x.lo)*(yr[1]-yr[0]) >= (y.hi-y.lo)*(xr[1]-xr[0]):
                mid = (x.lo+x.hi)//2
                assert x.lo < mid < x.hi
                stack.extend([(I(mid, x.hi), y, depth+1), (I(x.lo, mid), y, depth+1)])
            else:
                mid = (y.lo+y.hi)//2
                assert y.lo < mid < y.hi
                stack.extend([(x, I(mid, y.hi), depth+1), (x, I(y.lo, mid), depth+1)])
        now = time.monotonic()
        if now-last >= 10:
            print(f"{name}: visited={visited:,} leaves={sum(counts.values()):,} "
                  f"pending={len(stack)} covered={area/initial_area:.4%} "
                  f"elapsed={now-start:.1f}s", flush=True)
            last = now
    assert area == initial_area, "coverage accounting failed"
    assert sum(counts.values()) == split+len(roots)
    result = {
        "status": "PASS", "visited_boxes": visited, "split_boxes": split,
        "leaf_boxes": sum(counts.values()), "initial_boxes": len(roots),
        "max_depth": depth_max, "leaf_counts": dict(sorted(counts.items())),
        "coverage_exact": True, "leaf_sha256": digest.hexdigest(),
        "minimum_certified_margins": {
            k: {"exact": exact_str(v), "approx": v/Q}
            for k, v in sorted(margins.items())},
        "elapsed_seconds": round(time.monotonic()-start, 3),
    }
    print(name+": "+json.dumps(result, sort_keys=True), flush=True)
    return result


def self_test():
    """Check arithmetic enclosure against Fraction, including negative values."""
    samples = [I.rat(-7, 3), I.rat(-1, 10), I(-Q, Q//3),
               I.rat(0), I.rat(1, 7), I.rat(5, 3), I(Q//9, 2*Q)]
    def contains(box, f):
        assert Fraction(box.lo, Q) <= f <= Fraction(box.hi, Q)
    for x in samples:
        for a in (x.lo, (x.lo+x.hi)//2, x.hi):
            f = Fraction(a, Q)
            contains(x.sq(), f*f)
            for y in samples:
                for b in (y.lo, (y.lo+y.hi)//2, y.hi):
                    g = Fraction(b, Q)
                    contains(x+y, f+g)
                    contains(x-y, f-g)
                    contains(x*y, f*g)
                    if not y.lo <= 0 <= y.hi:
                        contains(x/y, f/g)
        if x.lo >= 0:
            r = x.sqrt()
            assert r.lo*r.lo <= x.lo*Q
            assert r.hi*r.hi >= x.hi*Q
    for n in (0, 1, 2, 3, 4, 1000000):
        x = I.rat(n)
        y = x.sqrt()
        assert y.lo*y.lo <= n*Q*Q <= y.hi*y.hi
    assert rho(I(-Q//2, 3*Q//5)).lo == 0
    assert rho(I(-Q//2, 3*Q//5)).hi == I.rat(1, 3).hi
    print("arithmetic_self_test: PASS", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "initialization", "drift"), default="all")
    parser.add_argument("--output", default="results.json")
    parser.add_argument("--max-boxes", type=int, default=5000000)
    parser.add_argument("--max-depth", type=int, default=70)
    args = parser.parse_args()
    self_test()
    results = {
        "status": "RUNNING", "arithmetic": "integer dyadic outward rounding",
        "bits": BITS, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": sys.version, "part": args.part,
        "parameters": {"c_star": "250001/1000000", "c_drift": ["1/4", "250001/1000000"],
                       "T": "53/25", "initialization_margin_target": "1/5000000",
                       "drift_increment": "1/1000", "drift_numerator_margin_target": "1/50000"},
    }
    path = Path(args.output)
    path.write_text(json.dumps(results, indent=2)+"\n")
    try:
        if args.part in ("all", "initialization"):
            d = CSTAR/(1+CSTAR)
            results["initialization"] = run_cover("initialization", initialization,
                (0, d.hi), (I.rat(-1, 2).lo, I.rat(3, 5).hi), args.max_boxes, args.max_depth)
            path.write_text(json.dumps(results, indent=2)+"\n")
        if args.part in ("all", "drift"):
            rng = (I.rat(-1, 2).lo, I.rat(3, 5).hi)
            results["drift"] = run_cover("drift", drift, rng, rng, args.max_boxes, args.max_depth)
        results["status"] = "PASS"
    except Exception as exc:
        results["status"] = "FAIL_OR_UNRESOLVED"
        results["error"] = str(exc)
        raise
    finally:
        path.write_text(json.dumps(results, indent=2)+"\n")
    print("ALL REQUESTED CHECKS: PASS", flush=True)


if __name__ == "__main__":
    main()
