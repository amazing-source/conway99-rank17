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
#include <time.h>

#define NV 49
#define CL 6
#ifndef DMIN
#define DMIN 2
#endif
#ifndef DMAX
#define DMAX 3
#endif
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
/* exact check of an integer candidate */
static int exact_neg(int n, const int64_t *M, const double *zd, int64_t *zout) {
    double mx = 0; for (int i = 0; i < n; i++) if (fabs(zd[i]) > mx) mx = fabs(zd[i]);
    if (mx == 0) return 0;
    for (double scale = 1e3; scale <= 1e7; scale *= 10) {
        int64_t z[25]; __int128 val = 0;
        for (int i = 0; i < n; i++) z[i] = (int64_t)llround(zd[i] / mx * scale);
        for (int i = 0; i < n; i++) { __int128 s = 0; for (int j = 0; j < n; j++) s += (__int128)M[i * n + j] * z[j]; val += (__int128)z[i] * s; }
        if (val < 0) { if (zout) memcpy(zout, z, n * sizeof(int64_t)); return 1; }
    }
    return 0;
}
/* LDL^T with diagonal pivoting: returns 1 with a candidate direction zd if a negative direction is found, 0 if PSD (tol) */
static int ldl_negdir(int n, const int64_t *M, double *zd) {
    double A[625]; int perm[25];
    for (int i = 0; i < n * n; i++) A[i] = (double)M[i];
    for (int i = 0; i < n; i++) perm[i] = i;
    const double tol = 1e-9;
    int i;
    double wv[25]; int fail = 0;
    for (i = 0; i < n; i++) {
        int p = i; for (int j = i + 1; j < n; j++) if (A[j * n + j] > A[p * n + p]) p = j;
        if (A[p * n + p] <= tol) {
            /* remaining Schur complement has diagonal <= tol */
            int jm = i; for (int j = i + 1; j < n; j++) if (A[j * n + j] < A[jm * n + jm]) jm = j;
            for (int t = 0; t < n; t++) wv[t] = 0;
            if (A[jm * n + jm] < -tol) { wv[jm] = 1; fail = 1; }
            else {
                double best = tol; int bj = -1, bk = -1;
                for (int j = i; j < n; j++) for (int k = j + 1; k < n; k++) if (fabs(A[j * n + k]) > best) { best = fabs(A[j * n + k]); bj = j; bk = k; }
                if (bj >= 0) { wv[bj] = 1; wv[bk] = (A[bj * n + bk] > 0) ? -1 : 1; fail = 1; }
            }
            break;
        }
        if (p != i) {
            for (int k = 0; k < n; k++) { double t = A[i * n + k]; A[i * n + k] = A[p * n + k]; A[p * n + k] = t; }
            for (int k = 0; k < n; k++) { double t = A[k * n + i]; A[k * n + i] = A[k * n + p]; A[k * n + p] = t; }
            int t = perm[i]; perm[i] = perm[p]; perm[p] = t;
        }
        double d = A[i * n + i];
        for (int j = i + 1; j < n; j++) {
            double l = A[j * n + i] / d;
            for (int k = i + 1; k < n; k++) A[j * n + k] -= l * A[i * n + k];
            A[j * n + i] = l;
        }
    }
    if (!fail) return 0;
    /* back-substitute x with L^T x = w (L unit lower in columns < i), then unpermute */
    double x[25];
    for (int t = n - 1; t >= 0; t--) {
        double s = wv[t];
        if (t < i) for (int j = t + 1; j < n; j++) s -= A[j * n + t] * x[j];
        x[t] = s;
    }
    for (int t = 0; t < n; t++) zd[perm[t]] = x[t];
    return 1;
}
static int neg_certificate(int n, const int64_t *M, int64_t *zout) {
    double zd[25];
    if (!ldl_negdir(n, M, zd)) return 0;
    if (exact_neg(n, M, zd, zout)) return 1;
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
static long long cnt_survivor = 0, cnt_nocert = 0, cnt_Qlevel_kill = 0;
static long long dnodes[6], dprune[6], d4_hsum = 0, d4_add = 0, d4_feas = 0, d4_ok = 0;

/* ---------------- M_U pieces (label level) ---------------- */
/* star label u = (j, s) <-> index u = 2j + (s<0); rep vertex starv(+1, j, s).  W label k <-> v_k. */
static int lab_rep(int u) { return (u < 12) ? starv(1, u / 2, (u % 2 == 0) ? 1 : -1) : vW(u - 12); }
static int Nlab(int u, int w) { int p = lab_rep(u), q = lab_rep(w); return hasE(p, sig(q)) - hasE(p, q); }
static void lab_root(int u, int r[7]) { root_vec(lab_rep(u), r); }
static int dot7(const int *a, const int *b) { int s = 0; for (int t = 0; t < 7; t++) s += a[t] * b[t]; return s; }

/* ---- bound tests.  Odd: M_U <= Mbound = [[K, M_SW],[M_WS, M_WW]] whenever K - M_SS >= 0; negative eigenvalue of
   Mbound restricted to S u {assigned W's} => M_U not PSD for every completion.  K0 = 20I - J (Q,P unknown),
   KQ = 15I + L - 4 Pi_Q - J (Q known, P unknown).  Even: form(f) <= F0(f) + 9*deg_unknown * sum_star f^2. ---- */
static int orbit_rep[NV], n_odd = 0, n_even = 0, odd_idx[NV], even_idx[NV];
static int Wmask = 31;                 /* bit k set: W label k has its star incidences fixed */
static int Qknown = 0;                 /* Q matching fixed */
static int Qp[NV], Tau[NV];
static int bound_neg_odd(void) {
    int R[17][7]; for (int u = 0; u < 17; u++) lab_root(u, R[u]);
    int lab[17], n = 0; for (int u = 0; u < 12; u++) lab[n++] = u;
    for (int k = 0; k < 5; k++) if (Wmask >> k & 1) lab[n++] = 12 + k;
    int64_t M[289];
    for (int ii = 0; ii < n; ii++) for (int jj = 0; jj < n; jj++) {
        int i = lab[ii], j = lab[jj];
        int64_t v;
        if (i < 12 && j < 12) {
            if (!Qknown) v = (i == j ? 20 : 0) - 1;
            else {
                int piQ = (Qp[lab_rep(i)] >= 0 && Qp[lab_rep(i)] == lab_rep(j)) ? 1 : 0;
                v = (i == j ? 15 : 0) + ((i ^ 1) == j ? 1 : 0) - 4 * piQ - 1;
            }
        } else v = 4 * Nlab(i, j) + (i == j ? 12 : 0) - dot7(R[i], R[j]);
        M[ii * n + jj] = v;
    }
    int r = neg_certificate(n, M, NULL);
    if (r < 0) { cnt_nocert++; return 0; }
    return r;
}
static int bound_neg_even(void) {
    /* vertices with known adjacency: x, inner, star, assigned W's (v_k, s v_k for k < Wassigned) */
    int use[NV]; int map[NV]; int n = 0;
    for (int v = 0; v < NV; v++) {
        use[v] = (v < 39) || (v < 44 ? (Wmask >> (v - 39) & 1) : (Wmask >> (v - 44) & 1));
    }
    for (int v = 0; v < NV; v++) map[v] = -1;
    int omap[25]; for (int i = 0; i < 25; i++) omap[i] = -1;
    for (int v = 0; v < NV; v++) if (use[v]) { int o = even_idx[v]; if (omap[o] < 0) omap[o] = n++; map[v] = omap[o]; }
    int64_t F[625]; memset(F, 0, sizeof(int64_t) * n * n);
    /* at the incidence level no star-star edge is present yet, at the Q level exactly the Q edges are present:
       the unknown star-star edges are the absent ones, degu per star vertex */
    for (int a = 0; a < NV; a++) if (use[a]) for (int b = 0; b < NV; b++) if (use[b]) {
        int A = hasE(a, b);
        F[map[a] * n + map[b]] += (a == b ? 27 : 0) - 9 * A + 1;
    }
    int degu = Qknown ? 1 : 2;
    for (int v = 15; v < 39; v++) F[map[v] * n + map[v]] += 9 * degu;   /* per vertex: -9 f^T A_unknown f <= 9 degu sum f^2 */
    int r = neg_certificate(n, F, NULL);
    if (r < 0) { cnt_nocert++; return 0; }
    return r;
}

/* ---------------- leaf tests (full Y' adjacency known) ---------------- */
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
static FILE *fdump = NULL; static long long ndump = 0, maxdump = 0;
static int leaf_tests(void) {
    static int64_t Pe[625], Po[625], Qe[625], Qo[625];
    if (fdump && ndump < maxdump) {
        ndump++;
        for (int k = 0; k < 5; k++) for (int t = 0; t < 6; t++) fprintf(fdump, "%d ", rho[k][t]);
        for (int v = 0; v < NV; v++) fprintf(fdump, "%llx ", (unsigned long long)adj[v]);
        for (int v = 15; v < 39; v++) fprintf(fdump, "%d ", hval[v]);
        fprintf(fdump, "\n"); fflush(fdump);
    }
    {   /* fast path: label-level M_U = 4 N_UU + 12 I - (R R^T)_UU, 17 x 17 */
        int R[17][7]; for (int u = 0; u < 17; u++) lab_root(u, R[u]);
        int64_t MU[289];
        for (int i = 0; i < 17; i++) for (int j = 0; j < 17; j++)
            MU[i * 17 + j] = (i == j ? 12 : 4 * Nlab(i, j)) - dot7(R[i], R[j]);
        int rr = neg_certificate(17, MU, NULL);
        if (rr > 0) { cnt_T1++; return 1; }
        if (rank_mod(17, MU) > 15) { cnt_T1r++; return 1; }
    }
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
    if (i == 12) {
        Qknown = 1;
        int kill = bound_neg_odd();
        if (kill) cnt_Qlevel_kill++; else qp_tau(0);
        Qknown = 0;
        return;
    }
    int y = 15 + i;
    if (Qp[y] >= 0) { qp_Q(i + 1); return; }
    for (int z = y + 1; z < 27; z++) {
        if (Qp[z] >= 0) continue;
        /* sign pre-constraints: need = +-2 forces h(Q(y)) = sign(need) */
        if (need[y] == 2 && hval[z] != 1) continue;
        if (need[y] == -2 && hval[z] != -1) continue;
        if (need[z] == 2 && hval[y] != 1) continue;
        if (need[z] == -2 && hval[y] != -1) continue;
        int r = add_checked(y, z);
        if (r == 0) continue;
        Qp[y] = z; Qp[z] = y; Qp[sig(y)] = sig(z); Qp[sig(z)] = sig(y);
        qp_Q(i + 1);
        Qp[y] = Qp[z] = Qp[sig(y)] = Qp[sig(z)] = -1;
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
static int ord[5];
/* look-ahead: every unassigned W must still admit 2+2 star neighbours with the right h-sum, each edge individually OK */
static int lookahead_ok(void) {
    for (int l = 0; l < 5; l++) {
        if (Wmask >> l & 1) continue;
        int v = vW(l), cnt[2][2] = {{0, 0}, {0, 0}};    /* [side][h=+1?] */
        for (int y = 15; y < 39; y++) {
            if (hasE(v, y)) { cnt[y >= 27][hval[y] > 0]++; continue; }
            int r = add_checked(v, y);
            if (r == 0) continue;
            if (r == 1) rem_pair(v, y);
            cnt[y >= 27][hval[y] > 0]++;
        }
        int npos = (sigma_k[l] + 4) / 2, found = 0;     /* number of +1 among the 4 */
        for (int a = 0; a <= 2 && !found; a++) {          /* a = #(+1) on side + */
            int b = npos - a; if (b < 0 || b > 2) continue;
            if (a <= cnt[0][1] && 2 - a <= cnt[0][0] && b <= cnt[1][1] && 2 - b <= cnt[1][0]) found = 1;
        }
        if (!found) return 0;
    }
    return 1;
}
static int cand[5][4356][4], ncand[5];
static long long nshards = 1, shard = 0, shard_ctr = 0;
static int count_only = 0;
static long long only_ww = -1, only_eps = -1, ww_local = 0, eps_local = 0;
static void build_cands(void) {
    for (int w = 0; w < 5; w++) {
        ncand[w] = 0;
        int v = vW(w);
        for (int i1 = 0; i1 < 12; i1++) for (int i2 = i1 + 1; i2 < 12; i2++)
        for (int i3 = 0; i3 < 12; i3++) for (int i4 = i3 + 1; i4 < 12; i4++) {
            int y[4] = {Op[i1], Op[i2], Om[i3], Om[i4]};
            if (hval[y[0]] + hval[y[1]] + hval[y[2]] + hval[y[3]] != sigma_k[w]) continue;
            int ok = 1, added[4], na = 0;
            for (int q = 0; q < 4 && ok; q++) { int r = add_checked(v, y[q]); if (r == 0) ok = 0; else added[na++] = (r == 1) ? y[q] : -1; }
            for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(v, added[q]);
            if (!ok) continue;
            memcpy(cand[w][ncand[w]++], y, sizeof y);
        }
    }
}
static long long knuth_walks = 0;
static double knuth_walk(int k, unsigned long long *rs) {
    /* returns estimate (product of branching counts) of wsum-ok leaves under this node along one random path */
    if (k >= 1 && k < 5 && !lookahead_ok()) return 0;
    if (k == 5) {
        for (int y = 15; y < 39; y++) { int t = wsum[y] + 4 * hval[y]; if (t < -2 || t > 2) return 0; }
        return 1;
    }
    int w = ord[k], v = vW(w), rem = 4 - k;
    static int ok_idx[5][4356]; int nok = 0;
    for (int c = 0; c < ncand[w]; c++) {
        int *y = cand[w][c]; int added[4], na = 0, ok = 1;
        for (int q = 0; q < 4; q++) { int r = add_checked(v, y[q]); if (r == 0) { ok = 0; break; } added[na++] = (r == 1) ? y[q] : -1; }
        if (ok) {
            for (int q = 0; q < 4; q++) { wsum[y[q]] += 2; wsum[sig(y[q])] -= 2; }
            int feas = 1;
            for (int s = 15; s < 39 && feas; s++) { int lo = -4 * hval[s] - 2, hi = -4 * hval[s] + 2; if (wsum[s] + 2 * rem < lo || wsum[s] - 2 * rem > hi) feas = 0; }
            for (int side = 0; side < 2 && feas; side++) {
                int pneed = 0, nneed = 0;
                for (int s = 15 + 12 * side; s < 27 + 12 * side; s++) { int lo = -4 * hval[s] - 2, hi = -4 * hval[s] + 2; if (wsum[s] < lo) pneed += (lo - wsum[s] + 1) / 2; if (wsum[s] > hi) nneed += (wsum[s] - hi + 1) / 2; }
                if (pneed > 2 * rem || nneed > 2 * rem) feas = 0;
            }
            if (feas) ok_idx[k][nok++] = c;
            for (int q = 0; q < 4; q++) { wsum[y[q]] -= 2; wsum[sig(y[q])] += 2; }
        }
        for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(v, added[q]);
    }
    if (nok == 0) return 0;
    *rs = *rs * 6364136223846793005ULL + 1442695040888963407ULL;
    int c = ok_idx[k][(*rs >> 33) % nok];
    int *y = cand[w][c]; int added[4], na = 0;
    for (int q = 0; q < 4; q++) { int r = add_checked(v, y[q]); added[na++] = (r == 1) ? y[q] : -1; }
    for (int q = 0; q < 4; q++) { wsum[y[q]] += 2; wsum[sig(y[q])] -= 2; }
    Wmask |= 1 << w;
    double sub = knuth_walk(k + 1, rs);
    Wmask &= ~(1 << w);
    for (int q = 0; q < 4; q++) { wsum[y[q]] -= 2; wsum[sig(y[q])] += 2; }
    for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(v, added[q]);
    return nok * sub;
}
static void inc_dfs(int k) {
    cnt_inc_nodes++;
    if (k >= 1 && k < 5 && !lookahead_ok()) { dprune[k]++; return; }
    if (k >= 4) {
        int kill = bound_neg_odd() || bound_neg_even();
        if (kill) { cnt_T0_prune++; dprune[k]++; return; }
    }
    if (k == 5) {
        cnt_inc_leaves++;
        for (int y = 15; y < 39; y++) { int t = wsum[y] + 4 * hval[y]; if (t < -2 || t > 2) return; }
        cnt_inc_wsum_ok++;
        run_QP();
        if ((cnt_inc_wsum_ok & 0x3FFF) == 0) { fprintf(stderr, "[%.0fs] inc_nodes %lld leaves %lld wsum_ok %lld T0 %lld Qlevel %lld qp %lld surv %lld | ww %lld eps %lld\n", (double)clock() / CLOCKS_PER_SEC, cnt_inc_nodes, cnt_inc_leaves, cnt_inc_wsum_ok, cnt_T0_prune, cnt_Qlevel_kill, cnt_qp_leaves, cnt_survivor, cnt_ww, cnt_eps); fflush(stderr); }
        return;
    }
    int w = ord[k], v = vW(w), rem = 4 - k;   /* W's after this one */
    for (int c = 0; c < ncand[w]; c++) {
        if (k == 0) { long long t = shard_ctr++; if (t % nshards != shard) continue; }
        int *y = cand[w][c];
        int added[4], na = 0, ok = 1;
        for (int q = 0; q < 4; q++) {
            int r = add_checked(v, y[q]);
            if (r == 0) { ok = 0; break; }
            added[na++] = (r == 1) ? y[q] : -1;
        }
        if (ok) {
            for (int q = 0; q < 4; q++) { wsum[y[q]] += 2; wsum[sig(y[q])] -= 2; }
            int feas = 1;
            for (int s = 15; s < 39 && feas; s++) {
                int lo = -4 * hval[s] - 2, hi = -4 * hval[s] + 2;
                if (wsum[s] + 2 * rem < lo || wsum[s] - 2 * rem > hi) feas = 0;
            }
            for (int side = 0; side < 2 && feas; side++) {
                int pneed = 0, nneed = 0;
                for (int s = 15 + 12 * side; s < 27 + 12 * side; s++) {
                    int lo = -4 * hval[s] - 2, hi = -4 * hval[s] + 2;
                    if (wsum[s] < lo) pneed += (lo - wsum[s] + 1) / 2;
                    if (wsum[s] > hi) nneed += (wsum[s] - hi + 1) / 2;
                }
                if (pneed > 2 * rem || nneed > 2 * rem) feas = 0;
            }
            if (feas) { Wmask |= 1 << w; inc_dfs(k + 1); Wmask &= ~(1 << w); }
            for (int q = 0; q < 4; q++) { wsum[y[q]] -= 2; wsum[sig(y[q])] += 2; }
        }
        for (int q = 0; q < na; q++) if (added[q] >= 0) rem_pair(v, added[q]);
    }
}

/* ---------------- W-W states and eps ---------------- */
static int PK[10], PL[10], ww_state[10];
/* ---- symmetry: stabilizer of the skeleton among signed permutations of lines 0..5 (x h-flip), acting on (ww state, eps) ---- */
static int ngrp = 0, gW[4096][5], gperm[4096][6], gflip[4096], PI_[5][5];
static int next_perm(int *a, int n) {
    int i = n - 2; while (i >= 0 && a[i] >= a[i + 1]) i--; if (i < 0) return 0;
    int j = n - 1; while (a[j] <= a[i]) j--; int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}
static void compute_stab(void) {
    ngrp = 0;
    int perm[6] = {0, 1, 2, 3, 4, 5};
    do {
        for (int sm = 0; sm < 64; sm++) {
            int img[5][6];
            for (int k = 0; k < 5; k++) for (int t = 0; t < 6; t++) img[k][perm[t]] = ((sm >> t) & 1 ? -1 : 1) * rho[k][t];
            for (int fl = 0; fl < 2; fl++) {
                int map[5], used = 0, ok = 1;
                for (int k = 0; k < 5 && ok; k++) {
                    int found = -1;
                    for (int l = 0; l < 5; l++) {
                        int eq = 1; for (int t = 0; t < 6; t++) if (img[k][t] != (fl ? -rho[l][t] : rho[l][t])) { eq = 0; break; }
                        if (eq) { found = l; break; }
                    }
                    if (found < 0 || (used >> found & 1)) ok = 0; else { map[k] = found; used |= 1 << found; }
                }
                if (ok && ngrp < 4096) {
                    memcpy(gW[ngrp], map, sizeof map); memcpy(gperm[ngrp], perm, sizeof perm); gflip[ngrp] = fl; ngrp++;
                }
            }
        }
    } while (next_perm(perm, 6));
}
/* canonical iff no group element (with or without c-flip) maps (state, eps) to a lexicographically smaller pair */
static int is_canonical(const int *state, const int *eps) {
    int base[16]; for (int p = 0; p < 10; p++) base[p] = state[p]; for (int t = 0; t < 6; t++) base[10 + t] = eps[t] + 1;
    for (int g = 0; g < ngrp; g++) for (int cf = 0; cf < 2; cf++) {
        int img[16];
        for (int p = 0; p < 10; p++) { int a = gW[g][PK[p]], b = gW[g][PL[p]]; img[PI_[a < b ? a : b][a < b ? b : a]] = state[p]; }
        for (int t = 0; t < 6; t++) { int e = eps[t] * (gflip[g] ? -1 : 1) * (cf ? -1 : 1); img[10 + gperm[g][t]] = e + 1; }
        for (int i = 0; i < 16; i++) { if (img[i] < base[i]) return 0; if (img[i] > base[i]) break; }
    }
    return 1;
}
static long long cnt_noncanon = 0;
static void eps_loop(void) {
    int net[6];
    for (int t = 0; t < 6; t++) { net[t] = 0; for (int k = 0; k < 5; k++) net[t] += rho[k][t]; }
    int sl[6], ns = 0;
    for (int t = 0; t < 6; t++) if (net[t] == 0) sl[ns++] = t;
    if (ns % 2) return;
    for (int mask = 0; mask < (1 << ns); mask++) {
        if (__builtin_popcount(mask) != ns / 2) continue;
        {   int epsv[6]; for (int t = 0; t < 6; t++) epsv[t] = 0;
            for (int q = 0; q < ns; q++) epsv[sl[q]] = (mask >> q & 1) ? 1 : -1;
            if (!is_canonical(ww_state, epsv)) { cnt_noncanon++; continue; } }
        cnt_eps++; eps_local++;
        if (only_eps >= 0 && eps_local - 1 != only_eps) continue;
        for (int t = 0; t < 6; t++) if (net[t] != 0) { a_lab[t][0] = -net[t]; a_lab[t][1] = net[t]; }
        for (int q = 0; q < ns; q++) { int e = (mask >> q & 1) ? 1 : -1; a_lab[sl[q]][0] = a_lab[sl[q]][1] = e; }
        for (int v = 15; v < 39; v++) { int al, j, be; st_root(v, &al, &j, &be); hval[v] = al * a_lab[j][(al * be == 1) ? 0 : 1]; }
        for (int y = 15; y < 39; y++) wsum[y] = 0;
        for (int q = 0; q < 5; q++) ord[q] = q;
        for (int x = 0; x < 5; x++) for (int yy = x + 1; yy < 5; yy++)
            if (abs(sigma_k[ord[yy]] + 4) > abs(sigma_k[ord[x]] + 4)) { int t = ord[x]; ord[x] = ord[yy]; ord[yy] = t; }
        Wmask = 0;
        build_cands();
        if (count_only) { long long tot = 1; for (int q = 0; q < 5; q++) tot *= ncand[q]; printf("  skel %lld wwl %lld epsl %lld ww %lld state %d%d%d%d%d%d%d%d%d%d eps %lld sig %d %d %d %d %d ncand %d %d %d %d %d (prod %.2e)\n", skel_index, ww_local - 1, eps_local - 1, cnt_ww, ww_state[0], ww_state[1], ww_state[2], ww_state[3], ww_state[4], ww_state[5], ww_state[6], ww_state[7], ww_state[8], ww_state[9], cnt_eps, sigma_k[0], sigma_k[1], sigma_k[2], sigma_k[3], sigma_k[4], ncand[0], ncand[1], ncand[2], ncand[3], ncand[4], (double)tot); }
        else if (knuth_walks > 0) {
            if (lookahead_ok()) {
                unsigned long long rs = 12345 + skel_index * 1000003ULL + eps_local * 7919ULL; double s = 0, s2 = 0;
                for (long long i = 0; i < knuth_walks; i++) { double e = knuth_walk(0, &rs); s += e; s2 += e * e; }
                double m = s / knuth_walks, sd = sqrt(fmax(s2 / knuth_walks - m * m, 0) / knuth_walks);
                printf("  KNUTH skel %lld wwl %lld epsl %lld est_leaves %.3e +- %.1e\n", skel_index, ww_local - 1, eps_local - 1, m, sd);
                fflush(stdout);
            }
        }
        else if (lookahead_ok()) inc_dfs(0);
        Wmask = 31;
    }
}
static void ww_dfs(int i, int *state) {
    if (i == 10) {
        for (int k = 0; k < 5; k++) {
            int al = 0, an = 0;
            for (int p = 0; p < 10; p++) if (PK[p] == k || PL[p] == k) { if (state[p] == 1) al++; if (state[p] == 2) an++; }
            int d = al - an; if (d < DMIN || d > DMAX) return;
            sigma_k[k] = -8 + 2 * d;
        }
        int nanti = 0, nD = 0; for (int p = 0; p < 10; p++) { if (state[p] == 2) nanti++; if (state[p] == 3) nD++; }
#ifdef NOANTI
        if (nanti > 0 || nD > 2) return;
#endif
        memcpy(ww_state, state, sizeof ww_state);
        cnt_ww++; ww_local++;
        if (only_ww >= 0 && ww_local - 1 != only_ww) return;
        eps_local = 0;
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
    if (argc > 6) { nshards = atoll(argv[5]); shard = atoll(argv[6]); }
    if (argc > 7) count_only = atoi(argv[7]);
    if (argc > 9) { fdump = fopen(argv[8], "w"); maxdump = atoll(argv[9]); }
    if (argc > 11) { only_ww = atoll(argv[10]); only_eps = atoll(argv[11]); }
    if (argc > 12) knuth_walks = atoll(argv[12]);
    int p = 0; for (int k = 0; k < 5; k++) for (int l = k + 1; l < 5; l++) { PK[p] = k; PL[p] = l; PI_[k][l] = p; p++; }
    for (int i = 0; i < 12; i++) { Op[i] = 15 + i; Om[i] = 27 + i; }
    setup_orbits();
    for (int s = 0; s < nsk; s++) {
        for (int k = 0; k < 5; k++) { for (int t = 0; t < 6; t++) if (fscanf(f, "%d", &rho[k][t]) != 1) return 1; rho[k][6] = 0; }
        if (s < first || s >= last) continue;
        skel_index = s;
        base_graph();
        compute_stab();
        /* the base graph must itself pass the partial check */
        int ok = 1; for (int v = 0; v < NV; v++) if (!vtx_ok(v)) ok = 0;
        if (!ok) { printf("skel %d: base graph inconsistent\n", s); continue; }
        int state[10];
        ww_local = 0;
        long long b_ww = cnt_ww, b_leaf = cnt_qp_leaves, b_surv = cnt_survivor, b_inc = cnt_inc_wsum_ok;
        ww_dfs(0, state);
        printf("skel %d: ww %lld incidence-leaves(wsum ok) %lld qp-leaves %lld survivors %lld | T0 kills %lld prunes %lld\n",
               s, cnt_ww - b_ww, cnt_inc_wsum_ok - b_inc, cnt_qp_leaves - b_leaf, cnt_survivor - b_surv, cnt_T0_kill, cnt_T0_prune);
        fflush(stdout);
    }
    printf("TOTAL ww %lld eps %lld inc_nodes %lld inc_leaves %lld wsum_ok %lld T0_prune %lld qp_leaves %lld | "
           "Qlevel %lld T1 %lld T1rank %lld T2 %lld T2rank %lld T3 %lld | survivors %lld nocert %lld\n",
           cnt_ww, cnt_eps, cnt_inc_nodes, cnt_inc_leaves, cnt_inc_wsum_ok, cnt_T0_prune, cnt_qp_leaves,
           cnt_Qlevel_kill, cnt_T1, cnt_T1r, cnt_T2, cnt_T2r, cnt_T3, cnt_survivor, cnt_nocert);
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

