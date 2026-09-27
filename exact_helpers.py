"""Exact, dependency-free helpers for the Conway involution frame search.

These verify certificates for supplied finite data. They do NOT certify shell
completeness, orbit completeness, or the graph-to-lattice reduction.

All coordinates in `binary_h_condition` must be IN AN INTEGRAL BASIS OF U,
not raw scaled Golay coordinates or arbitrary ambient Cartesian coordinates.
Run `python exact_helpers.py --self-test` for small arithmetic controls.
"""
from __future__ import annotations
from fractions import Fraction
from numbers import Integral
from typing import Sequence
import argparse


def _integer(x: object) -> int:
    if isinstance(x, Integral):
        return int(x)
    if isinstance(x, Fraction) and x.denominator == 1:
        return x.numerator
    raise TypeError(f"Expected an exact integer, got {x!r}")


def _matrix(a: Sequence[Sequence[int]], name: str, n: int | None = None) -> list[list[int]]:
    out = [[_integer(x) for x in row] for row in a]
    if n is None:
        if not out:
            raise ValueError(f"{name} must have at least one row")
        n = len(out[0])
    if any(len(row) != n for row in out):
        raise ValueError(f"{name} has inconsistent row lengths")
    return out


def box_farkas_gap(
    A: Sequence[Sequence[int]], b: Sequence[int], dual: Sequence[int],
    lower: Sequence[int] | None = None, upper: Sequence[int] | None = None,
    D: Sequence[Sequence[int]] | None = None,
    d: Sequence[int] | None = None,
    mu: Sequence[int] | None = None,
) -> int:
    """Return an EXACT gap; a strictly positive value certifies infeasibility.

    Problem: A x = b, D x <= d, lower <= x <= upper.
    Defaults: no D inequalities, and 0 <= x <= 1.
    Certificate: unrestricted integer dual, and nonnegative integer mu.

    Let c = A.T dual - D.T mu. Then every feasible x satisfies
      b.dual - d.mu <= c.x <= sum(max(lower[j]*c[j], upper[j]*c[j])).
    The returned gap is the leftmost quantity minus the rightmost one.

    Uses Python integers throughout, never fixed-width NumPy arithmetic.
    Rational data/certificates should first be cleared of denominators.
    """
    a = _matrix(A, "A")
    m, n = len(a), len(a[0])
    bv = [_integer(v) for v in b]
    lam = [_integer(v) for v in dual]
    if len(bv) != m or len(lam) != m:
        raise ValueError("b and dual must match the number of rows of A")
    lo = [0] * n if lower is None else [_integer(v) for v in lower]
    hi = [1] * n if upper is None else [_integer(v) for v in upper]
    if len(lo) != n or len(hi) != n or any(l > h for l, h in zip(lo, hi)):
        raise ValueError("Invalid bounds")
    c = [sum(a[i][j] * lam[i] for i in range(m)) for j in range(n)]
    target = sum(bv[i] * lam[i] for i in range(m))
    if D is not None:
        dd = _matrix(D, "D", n)
        if d is None or mu is None:
            raise ValueError("d and mu are required with D")
        dv, mv = [_integer(v) for v in d], [_integer(v) for v in mu]
        if len(dv) != len(dd) or len(mv) != len(dd) or any(v < 0 for v in mv):
            raise ValueError("Invalid inequality certificate")
        for i, row in enumerate(dd):
            target -= dv[i] * mv[i]
            for j in range(n):
                c[j] -= row[j] * mv[i]
    elif d is not None or mu is not None:
        raise ValueError("D is required when d or mu is supplied")
    support = sum(max(lo[j] * c[j], hi[j] * c[j]) for j in range(n))
    return target - support


def binary_rank(rows: Sequence[Sequence[int]]) -> int:
    a = _matrix(rows, "rows")
    pivots: dict[int, int] = {}
    for row in a:
        w = sum((x & 1) << j for j, x in enumerate(row))
        while w:
            p = w.bit_length() - 1
            if p in pivots:
                w ^= pivots[p]
            else:
                pivots[p] = w
                break
    return len(pivots)


def binary_h_condition(h_rows: Sequence[Sequence[int]]) -> bool:
    """Necessary binary-circuit condition at depths 1..7; not sufficient.

    Depth k<7: rank k. Depth 7: rank 6 AND total sum zero mod 2.
    Input rows are integral U-basis coefficients.
    """
    h = _matrix(h_rows, "h_rows")
    k = len(h)
    if not 1 <= k <= 7:
        raise ValueError("Expected 1 through 7 selected h-vectors")
    r = binary_rank(h)
    if k < 7:
        return r == k
    return r == 6 and all(sum(row[j] for row in h) % 2 == 0 for j in range(len(h[0])))


def short_dual_witness(
    G: Sequence[Sequence[int]], allowed_u: Sequence[Sequence[int]],
    y1: Sequence[int], y2: Sequence[int],
    numerator: Sequence[int], denominator: int,
) -> dict[str, object]:
    """Check w != 0, w perp y1,y2, norm(w)<8/7, (w,u) integral.

    G is a positive-definite integral Gram matrix in a fixed U-basis.
    w = numerator/denominator in that basis; allowed_u must include every u
    permitted by any completion of the search node. A valid witness excludes
    that node, USING the established min(N*) >= 8/7 theorem.

    The routine does not prove G positive definite or candidate completeness.
    """
    g = _matrix(G, "G")
    n = len(g)
    if any(len(row) != n for row in g):
        raise ValueError("G must be square")
    if any(g[i][j] != g[j][i] for i in range(n) for j in range(n)):
        raise ValueError("G must be symmetric")
    den = _integer(denominator)
    if den <= 0:
        raise ValueError("denominator must be positive")
    w = [_integer(x) for x in numerator]
    ys = [[_integer(x) for x in y1], [_integer(x) for x in y2]]
    us = _matrix(allowed_u, "allowed_u", n)
    if len(w) != n or any(len(y) != n for y in ys):
        raise ValueError("Coordinate dimension mismatch")
    gw = [sum(g[i][j] * w[j] for j in range(n)) for i in range(n)]
    norm_num = sum(w[i] * gw[i] for i in range(n))
    orthogonal = all(sum(y[i] * gw[i] for i in range(n)) == 0 for y in ys)
    integral = all(sum(u[i] * gw[i] for i in range(n)) % den == 0 for u in us)
    positive = norm_num > 0
    short = 7 * norm_num < 8 * den * den
    return {"valid": positive and short and orthogonal and integral,
            "norm": str(Fraction(norm_num, den * den)),
            "orthogonal_to_y": orthogonal, "integral_on_candidates": integral}


def self_test() -> None:
    # A feasible box has no certificate with positive gap.
    assert box_farkas_gap([[1, 1]], [1], [1]) == -1
    # Sum of two [0,1] variables cannot be 3.
    assert box_farkas_gap([[1, 1]], [3], [1]) == 1
    # Node bounds x0=x1=0 contradict x0+x1=1.
    assert box_farkas_gap([[1, 1]], [1], [1], upper=[0, 0]) == 1
    # x=1/2 together with x<=0.
    assert box_farkas_gap([[2]], [1], [1], D=[[1]], d=[0], mu=[2]) == 1
    # Explicit overflow control, well beyond int64.
    big = 2 ** 100
    assert box_farkas_gap([[big]], [big + 1], [1]) == 1
    e = [[int(i == j) for j in range(6)] for i in range(6)]
    assert binary_h_condition(e)
    assert binary_h_condition(e + [[1] * 6])
    assert not binary_h_condition([e[0], e[0]])
    # Arithmetic-only toy witness; not a Conway frame or a lattice control.
    test = short_dual_witness([[1,0,0],[0,1,0],[0,0,1]], [[2,0,0]],
                             [0,1,0], [0,0,1], [1,0,0], 2)
    assert test["valid"] and test["norm"] == "1/4"
    bad = short_dual_witness([[1,0,0],[0,1,0],[0,0,1]], [[1,0,0]],
                            [0,1,0], [0,0,1], [1,0,0], 2)
    assert not bad["valid"]
    # Exact finite-field controls for the even-weight model E_7.
    E = [r for r in range(128) if r.bit_count() % 2 == 0]
    q = lambda r: (r.bit_count() // 2) % 2
    assert len(E) == 64
    assert sum(q(r) == 0 for r in E) == 36  # plus type in dimension 6
    assert [r for r in E if all((r&s).bit_count()%2 == 0 for s in E)] == [0]
    roots = [(1<<i) | (1<<j) for i in range(7) for j in range(i+1,7)]
    for p in range(128):
        k = min(p.bit_count(), 7-p.bit_count())
        assert sum((p&r).bit_count()%2 for r in roots) == k*(7-k)
    # Four signed-permutation representatives for the rootless U, V=Z^17 subcase.
    for k in (0, 2, 4, 6):
        y1 = [1]*7 + [0]*10
        y2 = [1]*(k//2) + [-1]*(k//2) + [0]*(7-k) + [1]*(7-k) + [0]*(3+k)
        assert len(y2) == 17
        assert sum(a*a for a in y2) == 7 and sum(a*b for a,b in zip(y1,y2)) == 0
    # The other four Z17 pair types allowed by min(U)>=2: A-B and B-B.
    for overlap in (0, 2):
        for first_type in ("A", "B"):
            y1 = [0] * 17
            y2 = [0] * 17
            if first_type == "A":
                y1[:7] = [1] * 7
                y2[7] = 2
                common = list(range(overlap))
                external = list(range(8, 8 + 3 - overlap))
            else:
                y1[0] = 2
                y1[1:4] = [1] * 3
                y2[4] = 2
                common = list(range(1, 1 + overlap))
                external = list(range(5, 5 + 3 - overlap))
            for j, pos in enumerate(common):
                y2[pos] = 1 if j < overlap // 2 else -1
            for pos in external:
                y2[pos] = 1
            assert sum(a*a for a in y1) == sum(b*b for b in y2) == 7
            assert sum(a*b for a,b in zip(y1,y2)) == 0
            assert all(a*a+b*b <= 4 for a,b in zip(y1,y2))
    print("PASS: exact box certificates, inequality certificates, overflow control, "
          "binary circuit, dual witness, E7 parity model, eight Z17 y-pair types (four rootless)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        parser.print_help()
