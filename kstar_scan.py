"""Whole-pair exclusion scan via Sol's theorem (28/09): min(K*) >= 5/7 for every graph-origin K.
For each y-pair: H = 7 G^-1 (integral, since the discriminant group is 7-elementary). An integer z != 0 with z^T H z <= 4
gives kappa = G^-1 z in K* with 7|kappa|^2 = z^T H z <= 4 < 5, which excludes the whole pair. The certificate is b = H z with
G b = 7 z and b^T z <= 4, checked in exact integers.
Enumeration: LLL on H (floating point, only to pick a good basis), then Fincke-Pohst with bound 4 + eps. Every hit is
re-checked exactly. The minimum 7|kappa|^2 and the counts of vectors with 7|kappa|^2 = 1..6 are recorded.
Usage: python kstar_scan.py OUT.json FILE [FILE ...]"""
import sys, json, math
from fractions import Fraction
import numpy as np


def inv_times7(G):
    """exact 7*G^-1 via fraction Gauss-Jordan; asserts integrality"""
    n = len(G)
    A = [[Fraction(G[i][j]) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [x / piv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    H = [[7 * A[i][n + j] for j in range(n)] for i in range(n)]
    for row in H:
        for x in row:
            assert x.denominator == 1, 'H not integral'
    return [[int(x) for x in row] for row in H]


def lll(B, G, delta=0.99):
    """LLL on basis rows B (integer) for the quadratic form G; returns the reduced basis (integer rows)"""
    B = [list(r) for r in B]
    n = len(B)
    Gm = np.array(G, dtype=float)

    def ip(u, v):
        return float(np.array(u, dtype=float) @ Gm @ np.array(v, dtype=float))
    k = 1
    while k < n:
        # Gram-Schmidt (recomputed; n = 15 is small)
        Bs, mu, bb = [], [[0.0] * n for _ in range(n)], []
        for i in range(n):
            v = np.array(B[i], dtype=float)
            vi = v.copy()
            for j in range(i):
                mu[i][j] = float(v @ Gm @ Bs[j]) / bb[j]
                vi = vi - mu[i][j] * Bs[j]
            Bs.append(vi); bb.append(float(vi @ Gm @ vi))
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [x - q * y for x, y in zip(B[k], B[j])]
                for l in range(j + 1):
                    mu[k][l] -= q * (mu[j][l] if l < j else 1.0)
        if bb[k] >= (delta - mu[k][k - 1] ** 2) * bb[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            k = max(k - 1, 1)
    return B


def short_vectors(H, bound):
    """all z != 0 (up to sign) with z^T H z <= bound, exact; returns list of (value, z)"""
    n = len(H)
    B = lll([[int(i == j) for j in range(n)] for i in range(n)], H)
    Hr = [[sum(B[i][a] * H[a][b] * B[j][b] for a in range(n) for b in range(n)) for j in range(n)] for i in range(n)]
    M = np.array(Hr, dtype=float)
    L = np.linalg.cholesky(M)          # M = L L^T
    R = L.T                             # upper triangular
    out = []
    x = [0] * n
    eps = 1e-7

    def rec(i, rem):
        # coordinate i, given x[i+1..n-1]; the contribution of row i is (R[i,i] x_i + sum_{j>i} R[i,j] x_j)^2
        c = sum(R[i][j] * x[j] for j in range(i + 1, n))
        r = R[i][i]
        span = math.sqrt(max(rem, 0.0)) / r
        lo = math.ceil(-c / r - span - eps); hi = math.floor(-c / r + span + eps)
        for v in range(lo, hi + 1):
            x[i] = v
            t = (r * v + c) ** 2
            if t > rem + eps:
                continue
            if i == 0:
                if any(x):
                    zr = x[:]
                    z = [sum(zr[k] * B[k][a] for k in range(n)) for a in range(n)]
                    val = sum(z[a] * H[a][b] * z[b] for a in range(n) for b in range(n))
                    if val <= bound:
                        out.append((val, z))
            else:
                rec(i - 1, rem - t)
        x[i] = 0
    rec(n - 1, bound + 1e-6)
    # keep one of each +-z
    seen, res = set(), []
    for val, z in out:
        key = tuple(z) if tuple(z) > tuple(-a for a in z) else tuple(-a for a in z)
        if key in seen:
            continue
        seen.add(key); res.append((val, list(key)))
    return res


def main():
    outp = sys.argv[1]
    results = []
    for f in sys.argv[2:]:
        d = json.load(open(f, encoding='utf-8'))
        tag = f.split('K_')[-1].split('.json')[0]
        for rec_ in d['K']:
            G = rec_['K_gram']
            H = inv_times7(G)
            sv = short_vectors(H, 6)
            vals = sorted(v for v, _ in sv)
            mn = vals[0] if vals else None
            cnt = {k: sum(1 for v in vals if v == k) for k in range(1, 7)}
            wit = None
            if mn is not None and mn <= 4:
                z = min(sv)[1]
                n = len(G)
                b = [sum(H[i][j] * z[j] for j in range(n)) for i in range(n)]
                Gb = [sum(G[i][j] * b[j] for j in range(n)) for i in range(n)]
                assert Gb == [7 * zi for zi in z]
                assert sum(bi * zi for bi, zi in zip(b, z)) == min(sv)[0]
                wit = {'z': z, 'b': b, 'n': min(sv)[0]}
            results.append({'V': tag, 'id': rec_['id'], 'min7': mn, 'counts7': cnt, 'excluded': bool(mn is not None and mn <= 4),
                            'witness': wit})
        ex = sum(1 for r in results if r['V'] == tag and r['excluded'])
        tot = sum(1 for r in results if r['V'] == tag)
        print(f'{tag}: {ex}/{tot} excluded by min K* < 5/7', flush=True)
    json.dump(results, open(outp, 'w'), indent=0)
    print('total excluded', sum(r['excluded'] for r in results), 'of', len(results))


if __name__ == '__main__':
    main()
