/* r8k1.c -- exhaustive enumeration of the local configurations of R8, k = 1, with exact spectral kill tests.
 *
 * Setting (structures, 28/09): h is an integer sigma-odd (-4)-eigenvector of Gamma, zero on x u N(x), with
 *   h = +-1 on the 24 roots of star(c)  and  h = +-2 on the 10 roots of five labels W_1..W_5 off star(c).
 * Y' = {x} u N(x) u star(c) u W2 (49 vertices).  Everything on Y' that is forced by the roots is fixed; the unknowns are
 *   (1) the W-W adjacency, (2) which star roots each v_k sees (2 in O_c+, 2 in O_c-), (3) the matchings Q (inside O_c+-)
 *   and P (O_c+ <-> O_c-) of star(c).
 * Necessary conditions enforced: eigen-equation at every root of Y', lambda/mu upper bounds on all pairs of Y',
 *   p !~ sigma p, degrees <= 14.
 * Kill tests (all necessary for a genuine graph, since P_-4 = (27I - 9A + J)/63 and P_3 = (44I + 11A - 2J)/77 are
 *   orthogonal projectors of rank 44 and 54):
 *   T0 (crude, independent of Q,P): Sigma = M_WW - M_WS M_SS^+ M_SW must be PSD, and M_SS <= 20 on 1^perp, so
 *       20 M_WW - M_WS M_SW must be PSD.  [M_U = 4N_UU + 12I - (RR^T)_UU = 28 (Pi_4)_UU.]
 *   T1: odd part of 27I - 9A_Y' + J PSD, rank <= 22.     T2: even part PSD, rank <= 22.
 *   T3: odd and even parts of 44I + 11A_Y' - 2J PSD.
 * Every kill by a PSD test is certified by an integer vector z with z^T M z < 0 computed exactly in int64;
 * every kill by a rank test is certified by the rank modulo p = 2^31-1 (rank_p <= rank_Q).
 * Usage: r8k1 skeletons.txt first last   (skeleton index range, inclusive-exclusive)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define NV 49
#define CL 6
typedef uint64_t u64;

static int rho[5][7];
static int hval[NV];
static u64 adj[NV];
static int a_lab[6][2];          /* a on star label (j, s): index s: 0 -> s=+1, 1 -> s=-1 */

static inline int inner(int t, int s) { return 1 + 2 * t + (s == 1 ? 0 : 1); }
static inline int starv(int al, int j, int be) { return 15 + (al == 1 ? 0 : 12) + 2 * j + (be == 1 ? 0 : 1); }
static inline int vW(int k) { return 39 + k; }
static inline int sW(int k) { return 44 + k; }
static inline int sig(int v) {
    if (v == 0) return 0;
    if (v < 15) return ((v - 1) % 2 == 0) ? v + 1 : v - 1;
    if (v < 39) { int w = v - 15, al = w / 12, r = w % 12; return 15 + (1 - al) * 12 + (r ^ 1); }
    return (v < 44) ? v + 5 : v - 5;
}
static inline void st_root(int v, int *al, int *j, int *be) {
    int w = v - 15; *al = (w < 12) ? 1 : -1; w %= 12; *j = w / 2; *be = (w % 2 == 0) ? 1 : -1;
}
static void root_vec(int v, int r[7]) {
    memset(r, 0, 7 * sizeof(int));
    if (v >= 15 && v < 39) { int al, j, be; st_root(v, &al, &j, &be); r[CL] = al; r[j] = be; }
    else if (v >= 39 && v < 44) { for (int t = 0; t < 7; t++) r[t] = rho[v - 39][t]; }
    else if (v >= 44) { for (int t = 0; t < 7; t++) r[t] = -rho[v - 44][t]; }
}
static inline void addE(int a, int b) { adj[a] |= 1ULL << b; adj[b] |= 1ULL << a; }
static inline void remE(int a, int b) { adj[a] &= ~(1ULL << b); adj[b] &= ~(1ULL << a); }
static inline int hasE(int a, int b) { return (adj[a] >> b) & 1; }
static inline int pairok(int a, int b) {
    int c = __builtin_popcountll(adj[a] & adj[b]);
    return c <= (hasE(a, b) ? 1 : 2);
}
static int vtx_ok(int u) {
    if (__builtin_popcountll(adj[u]) > 14) return 0;
    for (int w = 0; w < NV; w++) if (w != u && !pairok(u, w)) return 0;
    return 1;
}
/* add edge a-b and its sigma image; keep it iff the partial lambda/mu check passes */
static int add_checked(int a, int b) {
    int sa = sig(a), sb = sig(b);
    if (a == sb) return 0;
    if (hasE(a, b)) return -1;               /* already there: nothing added */
    addE(a, b); addE(sa, sb);
    if (vtx_ok(a) && vtx_ok(b) && vtx_ok(sa) && vtx_ok(sb)) return 1;
    remE(a, b); remE(sa, sb); return 0;
}
static void rem_pair(int a, int b) { remE(a, b); remE(sig(a), sig(b)); }

/* ---------------- linear algebra ---------------- */
/* Jacobi eigenvalues of symmetric n x n (n <= 25); returns eigenvalues in w, eigenvectors in columns of V */
static void jacobi(int n, double *A, double *w, double *V) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) V[i * n + j] = (i == j);
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
        if (off < 1e-22) break;
        for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) {
            double apq = A[p * n + q];
            if (fabs(apq) < 1e-15) continue;
            double theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
            double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1));
            double c = 1 / sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < n; k++) {
                double akp = A[k * n + p], akq = A[k * n + q];
                A[k * n + p] = c * akp - s * akq; A[k * n + q] = s * akp + c * akq;
            }
            for (int k = 0; k < n; k++) {
                double apk = A[p * n + k], aqk = A[q * n + k];
                A[p * n + k] = c * apk - s * aqk; A[q * n + k] = s * apk + c * aqk;
            }
            for (int k = 0; k < n; k++) {
                double vkp = V[k * n + p], vkq = V[k * n + q];
                V[k * n + p] = c * vkp - s * vkq; V[k * n + q] = s * vkp + c * vkq;
            }
        }
    }
    for (int i = 0; i < n; i++) w[i] = A[i * n + i];
}
/* exact certificate: find integer z with z^T M z < 0 from the eigenvector of the smallest eigenvalue */
static int neg_certificate(int n, const int64_t *M, int64_t *zout) {
    double A[625], w[25], V[625];
    for (int i = 0; i < n * n; i++) A[i] = (double)M[i];
    jacobi(n, A, w, V);
    int im = 0; for (int i = 1; i < n; i++) if (w[i] < w[im]) im = i;
    if (w[im] > -1e-9) return 0;
    double mx = 0; for (int i = 0; i < n; i++) if (fabs(V[i * n + im]) > mx) mx = fabs(V[i * n + im]);
    for (double scale = 1e2; scale <= 1e6; scale *= 10) {
        int64_t z[25], val = 0;
        for (int i = 0; i < n; i++) z[i] = (int64_t)llround(V[i * n + im] / mx * scale);
        for (int i = 0; i < n; i++) { int64_t s = 0; for (int j = 0; j < n; j++) s += M[i * n + j] * z[j]; val += z[i] * s; }
        if (val < 0) { if (zout) memcpy(zout, z, n * sizeof(int64_t)); return 1; }
    }
    return -1;   /* numerically negative but no exact certificate found (should not happen) */
}
static const int64_t PMOD = 2147483647LL;
static int64_t pw(int64_t b, int64_t e) { int64_t r = 1; b %= PMOD; if (b < 0) b += PMOD; while (e) { if (e & 1) r = r * b % PMOD; b = b * b % PMOD; e >>= 1; } return r; }
static int rank_mod(int n, const int64_t *M) {
    int64_t A[625]; for (int i = 0; i < n * n; i++) { A[i] = M[i] % PMOD; if (A[i] < 0) A[i] += PMOD; }
    int r = 0;
    for (int c = 0; c < n && r < n; c++) {
        int piv = -1; for (int i = r; i < n; i++) if (A[i * n + c]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (int k = 0; k < n; k++) { int64_t t = A[r * n + k]; A[r * n + k] = A[piv * n + k]; A[piv * n + k] = t; }
        int64_t inv = pw(A[r * n + c], PMOD - 2);
        for (int k = 0; k < n; k++) A[r * n + k] = A[r * n + k] * inv % PMOD;
        for (int i = 0; i < n; i++) if (i != r && A[i * n + c]) {
            int64_t f = A[i * n + c];
            for (int k = 0; k < n; k++) { A[i * n + k] = (A[i * n + k] - f * A[r * n + k]) % PMOD; if (A[i * n + k] < 0) A[i * n + k] += PMOD; }
        }
        r++;
    }
    return r;
}

/* ---------------- counters ---------------- */
static long long cnt_ww = 0, cnt_eps = 0, cnt_inc_nodes = 0, cnt_inc_leaves = 0, cnt_inc_wsum_ok = 0;
static long long cnt_T0_kill = 0, cnt_T0_prune = 0, cnt_qp_leaves = 0, cnt_T1 = 0, cnt_T1r = 0, cnt_T2 = 0, cnt_T2r = 0, cnt_T3 = 0;
static long long cnt_survivor = 0, cnt_nocert = 0;

/* ---------------- M_U pieces (label level) ---------------- */
/* star label u = (j, s) <-> index u = 2j + (s<0); rep vertex starv(+1, j, s).  W label k <-> v_k. */
static int lab_rep(int u) { return (u < 12) ? starv(1, u / 2, (u % 2 == 0) ? 1 : -1) : vW(u - 12); }
static int Nlab(int u, int w) { int p = lab_rep(u), q = lab_rep(w); return hasE(p, sig(q)) - hasE(p, q); }
static void lab_root(int u, int r[7]) { root_vec(lab_rep(u), r); }
static int dot7(const int *a, const int *b) { int s = 0; for (int t = 0; t < 7; t++) s += a[t] * b[t]; return s; }

/* T0 on the first kk W labels: 20 M_WW - M_WS M_SW (leading block), needs the star-W incidences of those W's */
static int T0_block_neg(int kk) {
    int R[17][7]; for (int u = 0; u < 17; u++) lab_root(u, R[u]);
    int64_t T[25];
    int m[5][12];
    for (int k = 0; k < kk; k++) for (int u = 0; u < 12; u++) m[k][u] = 4 * Nlab(u, 12 + k) - dot7(R[u], R[12 + k]);
    for (int k = 0; k < kk; k++) for (int l = 0; l < kk; l++) {
        int mww = 4 * Nlab(12 + k, 12 + l) + (k == l ? 12 : 0) - dot7(R[12 + k], R[12 + l]);
        int64_t g = 0; for (int u = 0; u < 12; u++) g += (int64_t)m[k][u] * m[l][u];
        T[k * kk + l] = 20 * (int64_t)mww - g;
    }
    int r = neg_certificate(kk, T, NULL);
    if (r < 0) { cnt_nocert++; return 0; }
    return r;
}

/* ---------------- leaf tests (full Y' adjacency known) ---------------- */
static int orbit_rep[NV], n_odd = 0, n_even = 0, odd_idx[NV], even_idx[NV];
static void setup_orbits(void) {
    /* reps: x; inner (t,+); star with al=+1; v_k */
    n_odd = n_even = 0;
    for (int v = 0; v < NV; v++) { odd_idx[v] = even_idx[v] = -1; }
    for (int v = 0; v < NV; v++) {
        int s = sig(v);
        if (v == 0) { even_idx[0] = n_even++; continue; }
        if (v < s) { even_idx[v] = even_idx[s] = n_even++; odd_idx[v] = odd_idx[s] = n_odd++; orbit_rep[v] = 1; orbit_rep[s] = -1; }
    }
}
static void vertex_forms(int64_t *Pe, int64_t *Po, int64_t *Qe, int64_t *Qo) {
    /* P = 27I - 9A + J ; Q = 44I + 11A - 2J ; restricted to even / odd vectors in orbit coordinates */
    memset(Pe, 0, sizeof(int64_t) * n_even * n_even); memset(Po, 0, sizeof(int64_t) * n_odd * n_odd);
    memset(Qe, 0, sizeof(int64_t) * n_even * n_even); memset(Qo, 0, sizeof(int64_t) * n_odd * n_odd);
    for (int a = 0; a < NV; a++) for (int b = 0; b < NV; b++) {
        int A = hasE(a, b);
        int64_t p = (a == b ? 27 : 0) - 9 * A + 1, q = (a == b ? 44 : 0) + 11 * A - 2;
        int ea = even_idx[a], eb = even_idx[b];
        Pe[ea * n_even + eb] += p; Qe[ea * n_even + eb] += q;
        if (a != 0 && b != 0) {
            int oa = odd_idx[a], ob = odd_idx[b], s = orbit_rep[a] * orbit_rep[b];
            Po[oa * n_odd + ob] += s * p; Qo[oa * n_odd + ob] += s * q;
        }
    }
}
static int leaf_tests(void) {
    static int64_t Pe[625], Po[625], Qe[625], Qo[625];
    vertex_forms(Pe, Po, Qe, Qo);
    int r;
    r = neg_certificate(n_odd, Po, NULL); if (r > 0) { cnt_T1++; return 1; } if (r < 0) cnt_nocert++;
    if (rank_mod(n_odd, Po) > 22) { cnt_T1r++; return 1; }
    r = neg_certificate(n_even, Pe, NULL); if (r > 0) { cnt_T2++; return 1; } if (r < 0) cnt_nocert++;
    if (rank_mod(n_even, Pe) > 22) { cnt_T2r++; return 1; }
    r = neg_certificate(n_odd, Qo, NULL); if (r > 0) { cnt_T3++; return 1; }
    r = neg_certificate(n_even, Qe, NULL); if (r > 0) { cnt_T3++; return 1; }
    return 0;
}

/* ---------------- Q, P enumeration ---------------- */
static int Qp[NV], Tau[NV];     /* for star vertices of O_c+ (15..26): Q partner (in O_c+), tau partner (in O_c+) */
static int need[NV], wsum[NV];
static FILE *fsurv;
static long long skel_index;
static void dump_survivor(void) {
    cnt_survivor++;
    if (!fsurv) return;
    fprintf(fsurv, "SURVIVOR skel %lld rho", skel_index);
    for (int k = 0; k < 5; k++) { fprintf(fsurv, " ("); for (int t = 0; t < 6; t++) fprintf(fsurv, "%d%s", rho[k][t], t < 5 ? "," : ")"); }
    fprintf(fsurv, " adj");
    for (int v = 0; v < NV; v++) fprintf(fsurv, " %llx", (unsigned long long)adj[v]);
    fprintf(fsurv, "\n"); fflush(fsurv);
}
static void qp_tau(int i);
static void qp_Q(int i) {
    if (i == 12) { qp_tau(0); return; }
    int y = 15 + i;
    if (Qp[y] >= 0) { qp_Q(i + 1); return; }
    for (int z = y + 1; z < 27; z++) {
        if (Qp[z] >= 0) continue;
        int r = add_checked(y, z);
        if (r == 0) continue;
        Qp[y] = z; Qp[z] = y;
        qp_Q(i + 1);
        Qp[y] = Qp[z] = -1;
        if (r == 1) rem_pair(y, z);
    }
}
static void qp_tau(int i) {
    if (i == 12) {
        cnt_qp_leaves++;
        if (!leaf_tests()) dump_survivor();
        return;
    }
    int y = 15 + i;
    if (Tau[y] >= 0) { qp_tau(i + 1); return; }
    for (int z = y + 1; z < 27; z++) {
        if (Tau[z] >= 0) continue;
        /* sign rule: h(Q(y)) - h(tau(y)) = need(y), same for z */
        if (hval[Qp[y]] - hval[z] != need[y]) continue;
        if (hval[Qp[z]] - hval[y] != need[z]) continue;
        /* edges: y ~ sigma(z) and z ~ sigma(y) (and sigma images, added by add_checked) */
        int r1 = add_checked(y, sig(z));
        if (r1 == 0) continue;
        int r2 = add_checked(z, sig(y));
        if (r2 == 0) { if (r1 == 1) rem_pair(y, sig(z)); continue; }
        Tau[y] = z; Tau[z] = y;
        qp_tau(i + 1);
        Tau[y] = Tau[z] = -1;
        if (r2 == 1) rem_pair(z, sig(y));
        if (r1 == 1) rem_pair(y, sig(z));
    }
}
static void run_QP(void) {
    for (int v = 0; v < NV; v++) { Qp[v] = Tau[v] = -1; }
    for (int y = 15; y < 39; y++) need[y] = -4 * hval[y] - wsum[y];
    qp_Q(0);
}

/* ---------------- incidence enumeration ---------------- */
static int sigma_k[5];
static int Op[12], Om[12];
static void inc_dfs(int k) {
    cnt_inc_nodes++;
    if (k >= 2) {
        if (T0_block_neg(k)) { cnt_T0_prune++; return; }
    }
    if (k == 5) {
        cnt_inc_leaves++;
        for (int y = 15; y < 39; y++) { int t = wsum[y] + 4 * hval[y]; if (t < -2 || t > 2) return; }
        cnt_inc_wsum_ok++;
        run_QP();
        return;
    }
    int v = vW(k), rem = 4 - k;   /* W's after k */
    for (int i1 = 0; i1 < 12; i1++) for (int i2 = i1 + 1; i2 < 12; i2++)
    for (int i3 = 0; i3 < 12; i3++) for (int i4 = i3 + 1; i4 < 12; i4++) {
        int y[4] = {Op[i1], Op[i2], Om[i3], Om[i4]};
        if (hval[y[0]] + hval[y[1]] + hval[y[2]] + hval[y[3]] != sigma_k[k]) continue;
        int added[4], na = 0, ok = 1;
        for (int q = 0; q < 4; q++) {
            int r = add_checked(v, y[q]);
            if (r == 0) { ok = 0; break; }
            added[na++] = (r == 1) ? y[q] : -1;
        }
        if (ok) {
            for (int q = 0; q < 4; q++) { wsum[y[q]] += 2; wsum[sig(y[q])] -= 2; }
            /* feasibility of star W2-sums: remaining W's change each sum by at most 2 */
            int feas = 1;
            for (int s = 15; s < 39 && feas; s++) {
                int lo = -4 * hval[s] - 2, hi = -4 * hval[s] + 2;
                if (wsum[s] + 2 * rem < lo || wsum[s] - 2 * rem > hi) feas = 0;
            }
            if (feas) inc_dfs(k + 1);
            for (int q = 0; q < 4; q++) { wsum[y[q]] -= 2; wsum[sig(y[q])] += 2; }
        }
        for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(v, added[q]);
    }
}

/* ---------------- W-W states and eps ---------------- */
static int PK[10], PL[10];
static void eps_loop(void) {
    int net[6];
    for (int t = 0; t < 6; t++) { net[t] = 0; for (int k = 0; k < 5; k++) net[t] += rho[k][t]; }
    int sl[6], ns = 0;
    for (int t = 0; t < 6; t++) if (net[t] == 0) sl[ns++] = t;
    if (ns % 2) return;
    for (int mask = 0; mask < (1 << ns); mask++) {
        if (__builtin_popcount(mask) != ns / 2) continue;
        cnt_eps++;
        for (int t = 0; t < 6; t++) if (net[t] != 0) { a_lab[t][0] = -net[t]; a_lab[t][1] = net[t]; }
        for (int q = 0; q < ns; q++) { int e = (mask >> q & 1) ? 1 : -1; a_lab[sl[q]][0] = a_lab[sl[q]][1] = e; }
        for (int v = 15; v < 39; v++) { int al, j, be; st_root(v, &al, &j, &be); hval[v] = al * a_lab[j][(al * be == 1) ? 0 : 1]; }
        for (int y = 15; y < 39; y++) wsum[y] = 0;
        inc_dfs(0);
    }
}
static void ww_dfs(int i, int *state) {
    if (i == 10) {
        for (int k = 0; k < 5; k++) {
            int al = 0, an = 0;
            for (int p = 0; p < 10; p++) if (PK[p] == k || PL[p] == k) { if (state[p] == 1) al++; if (state[p] == 2) an++; }
            int d = al - an; if (d < 2 || d > 6) return;
            sigma_k[k] = -8 + 2 * d;
        }
        cnt_ww++;
        eps_loop();
        return;
    }
    int k = PK[i], l = PL[i];
    for (int s = 0; s < 4; s++) {
        int ea[3][2], ne = 0;
        if (s == 1 || s == 3) { ea[ne][0] = vW(k); ea[ne][1] = sW(l); ne++; ea[ne][0] = vW(l); ea[ne][1] = sW(k); ne++; }
        if (s == 2 || s == 3) { ea[ne][0] = vW(k); ea[ne][1] = vW(l); ne++; }
        int added[3], na = 0, ok = 1;
        for (int q = 0; q < ne; q++) {
            int r = add_checked(ea[q][0], ea[q][1]);
            if (r == 0) { ok = 0; break; }
            added[na++] = (r == 1) ? q : -1;
        }
        if (ok) {
            state[i] = s;
            int feas = 1;
            for (int m = 0; m < 5 && feas; m++) {
                int al = 0, an = 0, rem = 0;
                for (int p = 0; p < 10; p++) if (PK[p] == m || PL[p] == m) {
                    if (p <= i) { if (state[p] == 1) al++; if (state[p] == 2) an++; } else rem++;
                }
                if (al + rem - an < 2) feas = 0;
            }
            if (feas) ww_dfs(i + 1, state);
        }
        for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(ea[added[q]][0], ea[added[q]][1]);
    }
}
static void base_graph(void) {
    memset(adj, 0, sizeof adj);
    for (int t = 0; t < 7; t++) { addE(0, inner(t, 1)); addE(0, inner(t, -1)); addE(inner(t, 1), inner(t, -1)); }
    for (int v = 15; v < 49; v++) { int r[7]; root_vec(v, r); for (int t = 0; t < 7; t++) if (r[t]) addE(v, inner(t, r[t])); }
    for (int k = 0; k < 5; k++) { hval[vW(k)] = 2; hval[sW(k)] = -2; }
    for (int v = 0; v < 15; v++) hval[v] = 0;
}
int check_main(const char *fn);
int main(int argc, char **argv) {
    if (argc >= 3 && !strcmp(argv[1], "--check")) return check_main(argv[2]);
    if (argc < 4) { fprintf(stderr, "usage: r8k1 skeletons.txt first last [survivors.txt]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); int nsk; if (fscanf(f, "%d", &nsk) != 1) return 1;
    int first = atoi(argv[2]), last = atoi(argv[3]);
    fsurv = (argc > 4) ? fopen(argv[4], "w") : NULL;
    int p = 0; for (int k = 0; k < 5; k++) for (int l = k + 1; l < 5; l++) { PK[p] = k; PL[p] = l; p++; }
    for (int i = 0; i < 12; i++) { Op[i] = 15 + i; Om[i] = 27 + i; }
    setup_orbits();
    for (int s = 0; s < nsk; s++) {
        for (int k = 0; k < 5; k++) { for (int t = 0; t < 6; t++) if (fscanf(f, "%d", &rho[k][t]) != 1) return 1; rho[k][6] = 0; }
        if (s < first || s >= last) continue;
        skel_index = s;
        base_graph();
        /* the base graph must itself pass the partial check */
        int ok = 1; for (int v = 0; v < NV; v++) if (!vtx_ok(v)) ok = 0;
        if (!ok) { printf("skel %d: base graph inconsistent\n", s); continue; }
        int state[10];
        long long b_ww = cnt_ww, b_leaf = cnt_qp_leaves, b_surv = cnt_survivor, b_inc = cnt_inc_wsum_ok;
        ww_dfs(0, state);
        printf("skel %d: ww %lld incidence-leaves(wsum ok) %lld qp-leaves %lld survivors %lld | T0 kills %lld prunes %lld\n",
               s, cnt_ww - b_ww, cnt_inc_wsum_ok - b_inc, cnt_qp_leaves - b_leaf, cnt_survivor - b_surv, cnt_T0_kill, cnt_T0_prune);
        fflush(stdout);
    }
    printf("TOTAL ww %lld eps %lld inc_nodes %lld inc_leaves %lld wsum_ok %lld T0_prune %lld qp_leaves %lld | "
           "T1 %lld T1rank %lld T2 %lld T2rank %lld T3 %lld | survivors %lld nocert %lld\n",
           cnt_ww, cnt_eps, cnt_inc_nodes, cnt_inc_leaves, cnt_inc_wsum_ok, cnt_T0_prune, cnt_qp_leaves,
           cnt_T1, cnt_T1r, cnt_T2, cnt_T2r, cnt_T3, cnt_survivor, cnt_nocert);
    return 0;
}

/* ---- test harness: r8k1 --check file : each line = 30 ints (rho, lines 0..5) + 49 hex adjacency masks + 24 h values ---- */
int check_main(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) return 1;
    setup_orbits();
    int n = 0;
    while (1) {
        for (int k = 0; k < 5; k++) { for (int t = 0; t < 6; t++) if (fscanf(f, "%d", &rho[k][t]) != 1) goto done; rho[k][6] = 0; }
        for (int v = 0; v < NV; v++) { unsigned long long m; if (fscanf(f, "%llx", &m) != 1) goto done; adj[v] = m; }
        static int64_t Pe[625], Po[625], Qe[625], Qo[625];
        vertex_forms(Pe, Po, Qe, Qo);
        double A[625], w[25], V[625];
        double mo, me; int i;
        for (i = 0; i < n_odd * n_odd; i++) A[i] = (double)Po[i];
        jacobi(n_odd, A, w, V); mo = w[0]; for (i = 1; i < n_odd; i++) if (w[i] < mo) mo = w[i];
        for (i = 0; i < n_even * n_even; i++) A[i] = (double)Pe[i];
        jacobi(n_even, A, w, V); me = w[0]; for (i = 1; i < n_even; i++) if (w[i] < me) me = w[i];
        int ok = 1; for (int v = 0; v < NV; v++) if (!vtx_ok(v)) ok = 0;
        printf("sample %d: partial_ok %d odd_min %.4f rank_odd %d even_min %.4f rank_even %d\n", n, ok, mo, rank_mod(n_odd, Po), me, rank_mod(n_even, Pe));
        n++;
    }
done:
    return 0;
}
