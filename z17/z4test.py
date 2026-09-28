"""Exact mod-4 lift test for a Z^17 skeleton (structures, 28/09 ~03:15).

If the 42 frame codewords q_u = (r_u | w_u) are pairwise orthogonal mod 4, they generate (with 2 D^perp) a self-dual Z4 code Q with
residue D = graph of rho.  Self-dual codes with residue D <-> functions phi: D -> D* = F_2^24 / D^perp with
    phi(d + d') = phi(d) + phi(d') + (d * d')            (additivity; * = componentwise product)
    phi(d)(d') + phi(d')(d) = |d cap d'| / 2  (mod 2)    (self-orthogonality)
phi is fixed by M_c = phi(b_c) for the basis b_c = d_{0c} (the vertex-1 cells), c = 1..6: 21 free bits.
For a codeword q with residue d: sigma(q) = ((q - dtilde) mod 4 / 2) mod D^perp, read in the basis b (6 bits).
q in Q  <=>  sigma(q) = phi(d).  Cell cd needs an admissible pair (w+,w-) with sigma(q+) = sigma(q-) = phi(d_cd).
So the mod-4 relaxation (cross-cell rule mod 4 only; cell domains as built) is feasible  <=>  some M satisfies
phi(d_cd) in Sigma_cd for all 21 cells.  This is a necessary condition of the full problem (integer rules imply the mod-4 rule).
"""
import numpy as np
import stageB as B


def vec24(e_part, w):
    return np.concatenate([np.asarray(e_part, dtype=np.int64), np.asarray(w, dtype=np.int64)])


def z4_feasible(sk, return_info=False):
    pis = sk.pis
    # basis b_c (c = 1..6) of D as 24-bit int arrays: e-part eps_0 + eps_c, f-part pi_c (pi_0 = 0)
    b = []
    for c in range(1, 7):
        v = np.zeros(24, dtype=np.int64)
        v[0] = 1; v[c] = 1
        for i in range(17):
            if pis[c] >> i & 1:
                v[7 + i] = 1
        b.append(v)
    b = np.array(b)                                   # 6 x 24
    # J[c][e] = |b_c cap b_e| / 2 mod 2 ; kappa[c][d][e] = |b_c cap b_d cap b_e| mod 2
    J = ((b @ b.T) // 2) % 2
    Sig = {}
    for (c, d) in B.CELLS:
        P, Mn = sk.dom[(c, d)]
        if len(P) == 0:
            return (False, 'empty cell') if return_info else False
        ep = np.zeros(7, dtype=np.int64); ep[c] = 1; ep[d] = 1
        em = np.zeros(7, dtype=np.int64); em[c] = 1; em[d] = -1
        Qp = np.hstack([np.tile(ep, (len(P), 1)), P])
        Qm = np.hstack([np.tile(em, (len(Mn), 1)), Mn])
        def sig(Q):
            dt = Q % 2
            x = ((Q - dt) % 4) // 2
            s = (x @ b.T) % 2                          # n x 6
            return s @ (1 << np.arange(6))
        sp, sm = sig(Qp), sig(Qm)
        assert (sp == sm).all(), 'M1 should give equal signatures'
        Sig[(c, d)] = set(int(x) for x in sp)
    # kappa for cells (c,d), 1<=c<d<=6 (in basis indices c-1, d-1)
    kap = {}
    for c in range(1, 7):
        for d in range(c + 1, 7):
            val = 0
            for e in range(1, 7):
                if int((b[c - 1] * b[d - 1] * b[e - 1]).sum()) % 2:
                    val |= 1 << (e - 1)
            kap[(c, d)] = val
    # DFS over M_1..M_6 (6-bit rows, bit e-1 = phi(b_c)(b_e))
    M = [0] * 7
    sols = [0]

    def rec(c):
        if c == 7:
            sols[0] += 1
            return True
        for m in Sig[(0, c)]:
            ok = True
            for e in range(1, c):                      # self-orthogonality: M_c[e] + M_e[c] = J[c][e]
                if ((m >> (e - 1)) & 1) ^ ((M[e] >> (c - 1)) & 1) != J[c - 1][e - 1]:
                    ok = False; break
            if not ok:
                continue
            for e in range(1, c):                      # cells (e, c): phi = M_e + M_c + kappa
                if (M[e] ^ m ^ kap[(e, c)]) not in Sig[(e, c)]:
                    ok = False; break
            if not ok:
                continue
            M[c] = m
            if rec(c + 1):
                return True
        return False

    feas = rec(1)
    if return_info:
        return feas, {cd: len(Sig[cd]) for cd in B.CELLS}
    return feas


if __name__ == '__main__':
    import random, sys, time
    for p, fn in ((0, 'out/stageA_p0.raw'), (4, 'out/c_leaves_p4_m1.txt')):
        z = B.Z17(p)
        lines = open(fn).read().split()
        random.seed(21)
        samp = random.sample(lines, 40)
        n = 0; res = {}; t = time.time()
        for l in samp:
            cols, pis = B.parse_raw(l)
            sk = B.Skeleton(z, pis).build()
            if min(sk.sizes()) == 0:
                continue
            n += 1
            if n > 6:
                break
            f, info = z4_feasible(sk, True)
            res[f] = res.get(f, 0) + 1
        print('pair', p, 'z4-feasible counts:', res, '%.1fs' % (time.time() - t), flush=True)
