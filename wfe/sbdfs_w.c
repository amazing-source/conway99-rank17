/* sbdfs_w: sbdfs_g + witness budgets (structures, 28/09 night).  Extra optional input section after the 21 cells:
     int32 ncons; per constraint: int32 target; then per cell x (CELLS order) ndom[x] int8 contributions.
   A constraint says: sum over the 21 cells of contrib[cell][chosen entry] == target (witness Parseval sums a^2 = 4n and
   cross sums a(kappa)a(kappa') = 28(kappa,kappa')).  Pruning: running sum + suffix min/max (static over the domains, or
   dynamic over the live domains when argv[2] = 1).
   Original header:
   Stage-B exact DFS for a general rank-17 unimodular V (structures, 28/09).  Rules as in SPEC_ENGINE.md / sbdfs3.c:
   cells in lexicographic order (c<d, 0-based); caller-supplied domains of ordered signed pairs (w+, w-) already filtered by
   (w+,w-) in {0,+-4}, M1, refined L1 and sign normalisation; pairwise rule ((w_s^X, w_t^Y) + r_s.r_t)/4 in {0,+-1} with
   precomputed compatibility bitsets; star closure; star look-ahead with 1 or 2 open cells; PSD residual (Cholesky, 15 dims).
   Vectors are integer coordinates in a basis of V with Gram G; inner products are exact: (a,b) = a^T G b.
   Input (binary, little endian):
     int32 cap; int32 G[17*17]; double M[17*15]  (projection to an orthonormal basis of W: proj(x)_j = sum_i x_i M[i][j]);
     21 cells in order: int32 n; n records of 34 int8 (w+[17], w-[17]).
   Output: "VERDICT nodes d0..d21" (+ the 21 chosen records if FOUND).                                                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define NC 21
typedef uint64_t u64;
static int CC[NC], CD[NC], ndom[NC], nw[NC];
static signed char *dom[NC];
static int *gdom[NC];                  /* n*34 : G w+ and G w- */
static double *proj[NC];
static u64 *kp[NC], *km[NC];           /* 3 words per entry: packed 17 bytes of w+ + w- / w+ - w- */
static u64 *compat[NC][NC];
static u64 *lev[NC + 1][NC];
static int G[17][17];
static double Mm[17 * 15];
static long long nodes = 0, cap, depthc[NC + 1];
static int chosen[NC], rrt[NC][NC][2][2], lastk[7];
static int *htab; static int hcap = 0;   /* stores entry indices + 1 */
#define MAXCONS 2048
static int ncons = 0, ctarget[MAXCONS], dyn = 0;
static signed char *ccon[MAXCONS][NC];
static int csmin[MAXCONS][NC + 1], csmax[MAXCONS][NC + 1];
static int crun[NC + 1][MAXCONS];
static long long cons_pruned = 0, dyn_pruned = 0;

static int dyn_ok(int k) {                  /* live-domain bounds for the cells after k (lev[k + 1]) */
    for (int i = 0; i < ncons; i++) {
        int lo = crun[k + 1][i], hi = crun[k + 1][i];
        for (int Y = k + 1; Y < NC; Y++) {
            int mn = 1 << 30, mx = -(1 << 30);
            const u64 *b = lev[k + 1][Y];
            const signed char *cc = ccon[i][Y];
            for (int w = 0; w < nw[Y]; w++) {
                u64 x = b[w];
                while (x) {
                    int bit = __builtin_ctzll(x); x &= x - 1;
                    int v = cc[64 * w + bit];
                    if (v < mn) mn = v;
                    if (v > mx) mx = v;
                }
            }
            lo += mn; hi += mx;
        }
        if (lo > ctarget[i] || hi < ctarget[i]) return 0;
    }
    return 1;
}

static int pack(const int *v, u64 *out) {
    unsigned char b[24];
    memset(b, 0, sizeof b);
    for (int i = 0; i < 17; i++) {
        if (v[i] < -127 || v[i] > 127) return 0;
        b[i] = (unsigned char)(v[i] + 128);
    }
    memcpy(out, b, 24);
    return 1;
}

static u64 hash3(const u64 *k) {
    u64 h = k[0] * 0x9E3779B97F4A7C15ULL ^ (k[1] + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL ^ k[2] * 0x165667B19E3779F9ULL;
    return h ^ (h >> 29);
}

static int psd_ok(double R[15][15]) {
    double L[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j <= i; j++) {
            double s = R[i][j] + (i == j ? 1e-6 : 0.0);
            for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
            if (i == j) { if (s <= 0) return 0; L[i][i] = sqrt(s); }
            else L[i][j] = s / L[j][j];
        }
    return 1;
}

static void starsum(int v, int k, int *s) {
    memset(s, 0, 17 * sizeof(int));
    for (int x = 0; x <= k; x++) {
        const signed char *w = dom[x] + 34 * chosen[x];
        if (CC[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] + w[17 + i];
        else if (CD[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] - w[17 + i];
    }
}

static int anybit(const u64 *b, int n) { for (int i = 0; i < n; i++) if (b[i]) return 1; return 0; }

static int solve_rec(int k, double R[15][15]);

static int try_entry(int k, int t, double R[15][15]) {
    for (int i = 0; i < ncons; i++) {
        int v = crun[k][i] + ccon[i][k][t];
        if (v + csmin[i][k + 1] > ctarget[i] || v + csmax[i][k + 1] < ctarget[i]) { cons_pruned++; return 0; }
        crun[k + 1][i] = v;
    }
    const double *pv = proj[k] + 30 * t;
    double R2[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pv[i] * pv[j] - pv[15 + i] * pv[15 + j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    chosen[k] = t;
    int s[17];
    for (int v = 0; v < 7; v++)
        if (lastk[v] == k) { starsum(v, k, s); for (int i = 0; i < 17; i++) if (s[i]) return 0; }
    for (int Y = k + 1; Y < NC; Y++) {
        const u64 *cm = compat[k][Y] + (size_t)t * nw[Y];
        u64 *src = lev[k][Y], *dst = lev[k + 1][Y], acc = 0;
        for (int i = 0; i < nw[Y]; i++) { dst[i] = src[i] & cm[i]; acc |= dst[i]; }
        if (!acc) return 0;
    }
    for (int v = 0; v < 7; v++) {
        int op[3], cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CC[Y] == v || CD[Y] == v) { if (cnt < 3) op[cnt] = Y; cnt++; }
        if (cnt == 1) {
            int Y = op[0];
            starsum(v, k, s);
            u64 *b = lev[k + 1][Y];
            for (int w = 0; w < nw[Y]; w++) {
                u64 x = b[w];
                while (x) {
                    int bit = __builtin_ctzll(x); x &= x - 1;
                    const signed char *u = dom[Y] + 34 * (64 * w + bit);
                    for (int i = 0; i < 17; i++) {
                        int c = (CC[Y] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                        if (c != -s[i]) { b[w] &= ~(1ULL << bit); break; }
                    }
                }
            }
            if (!anybit(b, nw[Y])) return 0;
        } else if (cnt == 2) {
            starsum(v, k, s);
            for (int side = 0; side < 2; side++) {
                int A = op[side], Bc = op[1 - side];
                u64 *bB = lev[k + 1][Bc], *bA = lev[k + 1][A];
                u64 *kb = (CC[Bc] == v) ? kp[Bc] : km[Bc];
                int nB = 0;
                for (int w = 0; w < nw[Bc]; w++) nB += __builtin_popcountll(bB[w]);
                int hs = 1;
                while (hs < 2 * nB + 2) hs <<= 1;
                if (hs > hcap) { free(htab); htab = malloc(sizeof(int) * (size_t)hs); hcap = hs; }
                memset(htab, 0, sizeof(int) * (size_t)hs);
                for (int w = 0; w < nw[Bc]; w++) {
                    u64 x = bB[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        int q = 64 * w + bit;
                        int pos = (int)(hash3(kb + 3 * q) & (u64)(hs - 1));
                        while (htab[pos]) pos = (pos + 1) & (hs - 1);
                        htab[pos] = q + 1;
                    }
                }
                for (int w = 0; w < nw[A]; w++) {
                    u64 x = bA[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        const signed char *u = dom[A] + 34 * (64 * w + bit);
                        int need[17];
                        for (int i = 0; i < 17; i++) {
                            int c = (CC[A] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                            need[i] = -s[i] - c;
                        }
                        u64 key[3];
                        int found = 0;
                        if (pack(need, key)) {
                            int pos = (int)(hash3(key) & (u64)(hs - 1));
                            while (htab[pos]) {
                                const u64 *kk = kb + 3 * (htab[pos] - 1);
                                if (kk[0] == key[0] && kk[1] == key[1] && kk[2] == key[2]) { found = 1; break; }
                                pos = (pos + 1) & (hs - 1);
                            }
                        }
                        if (!found) bA[w] &= ~(1ULL << bit);
                    }
                }
                if (!anybit(bA, nw[A])) return 0;
            }
        }
    }
    if (dyn && ncons && !dyn_ok(k)) { dyn_pruned++; return 0; }
    return solve_rec(k + 1, R2);
}

static int solve_rec(int k, double R[15][15]) {
    nodes++;
    depthc[k]++;
    if (nodes > cap) return 2;
    if (k == NC) {
        for (int i = 0; i < 15; i++)
            for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0;
        return 1;
    }
    u64 *b = lev[k][k];
    for (int w = 0; w < nw[k]; w++) {
        u64 x = b[w];
        while (x) {
            int bit = __builtin_ctzll(x); x &= x - 1;
            int r = try_entry(k, 64 * w + bit, R);
            if (r) return r;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    if (!f) { printf("ERR open\n"); return 1; }
    int capi;
    if (fread(&capi, 4, 1, f) != 1) return 1;
    cap = capi;
    if (fread(G, sizeof(int), 17 * 17, f) != 17 * 17) return 1;
    if (fread(Mm, sizeof(double), 17 * 15, f) != 17 * 15) return 1;
    int k = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { CC[k] = c; CD[k] = d; k++; }
    for (int x = 0; x < NC; x++) {
        if (fread(&ndom[x], 4, 1, f) != 1) return 1;
        int n = ndom[x] ? ndom[x] : 1;
        nw[x] = (n + 63) / 64;
        dom[x] = malloc(34 * (size_t)n);
        if (ndom[x] && fread(dom[x], 1, 34 * (size_t)ndom[x], f) != 34 * (size_t)ndom[x]) return 1;
        gdom[x] = malloc(sizeof(int) * 34 * (size_t)n);
        proj[x] = malloc(sizeof(double) * 30 * (size_t)n);
        kp[x] = malloc(sizeof(u64) * 3 * (size_t)n); km[x] = malloc(sizeof(u64) * 3 * (size_t)n);
        for (int t = 0; t < ndom[x]; t++) {
            const signed char *w = dom[x] + 34 * t;
            int vp[17], vm[17];
            for (int i = 0; i < 17; i++) { vp[i] = w[i] + w[17 + i]; vm[i] = w[i] - w[17 + i]; }
            if (!pack(vp, kp[x] + 3 * t) || !pack(vm, km[x] + 3 * t)) { printf("ERR range\n"); return 1; }
            for (int s2 = 0; s2 < 2; s2++) {
                for (int i = 0; i < 17; i++) {
                    int a = 0;
                    for (int j = 0; j < 17; j++) a += G[i][j] * w[17 * s2 + j];
                    gdom[x][34 * t + 17 * s2 + i] = a;
                }
                for (int j = 0; j < 15; j++) {
                    double a = 0;
                    for (int i = 0; i < 17; i++) a += w[17 * s2 + i] * Mm[i * 15 + j];
                    proj[x][30 * t + 15 * s2 + j] = a;
                }
            }
        }
        for (int L = 0; L <= NC; L++) lev[L][x] = calloc((size_t)nw[x], sizeof(u64));
        for (int t = 0; t < ndom[x]; t++) lev[0][x][t >> 6] |= 1ULL << (t & 63);
    }
    if (fread(&ncons, 4, 1, f) != 1) ncons = 0;
    if (ncons > MAXCONS) { puts("ERR ncons"); return 1; }
    for (int i = 0; i < ncons; i++) {
        if (fread(&ctarget[i], 4, 1, f) != 1) { puts("ERR cons"); return 1; }
        for (int x = 0; x < NC; x++) {
            ccon[i][x] = malloc(ndom[x] ? (size_t)ndom[x] : 1);
            if (ndom[x] && fread(ccon[i][x], 1, (size_t)ndom[x], f) != (size_t)ndom[x]) { puts("ERR cons"); return 1; }
        }
        csmin[i][NC] = csmax[i][NC] = 0;
        for (int x = NC - 1; x >= 0; x--) {
            int mn = 1 << 30, mx = -(1 << 30);
            for (int t = 0; t < ndom[x]; t++) { int v = ccon[i][x][t]; if (v < mn) mn = v; if (v > mx) mx = v; }
            if (!ndom[x]) mn = mx = 0;
            csmin[i][x] = csmin[i][x + 1] + mn; csmax[i][x] = csmax[i][x + 1] + mx;
        }
    }
    if (argc > 2) dyn = atoi(argv[2]);
    fclose(f);
    for (int X = 0; X < NC; X++)
        for (int Y = 0; Y < NC; Y++)
            for (int s = 0; s < 2; s++)
                for (int t = 0; t < 2; t++) {
                    int r1[7] = {0}, r2[7] = {0}, v = 0;
                    r1[CC[X]] = 1; r1[CD[X]] = s ? -1 : 1;
                    r2[CC[Y]] = 1; r2[CD[Y]] = t ? -1 : 1;
                    for (int i = 0; i < 7; i++) v += r1[i] * r2[i];
                    rrt[X][Y][s][t] = v;
                }
    int empty = 0;
    for (int x = 0; x < NC; x++) if (!ndom[x]) empty = 1;
    if (!empty)
        for (int X = 0; X < NC; X++)
            for (int Y = X + 1; Y < NC; Y++) {
                compat[X][Y] = calloc((size_t)ndom[X] * nw[Y], sizeof(u64));
                for (int t = 0; t < ndom[X]; t++) {
                    const signed char *a = dom[X] + 34 * t;
                    u64 *row = compat[X][Y] + (size_t)t * nw[Y];
                    for (int q = 0; q < ndom[Y]; q++) {
                        const int *gb = gdom[Y] + 34 * q;
                        int ok = 1;
                        for (int s = 0; s < 2 && ok; s++)
                            for (int tt = 0; tt < 2 && ok; tt++) {
                                int g = rrt[X][Y][s][tt];
                                for (int i = 0; i < 17; i++) g += a[17 * s + i] * gb[17 * tt + i];
                                if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                            }
                        if (ok) row[q >> 6] |= 1ULL << (q & 63);
                    }
                }
            }
    for (int v = 0; v < 7; v++) lastk[v] = -1;
    for (int x = 0; x < NC; x++) { lastk[CC[x]] = x; lastk[CD[x]] = x; }
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    for (int i = 0; i < ncons; i++) {
        crun[0][i] = 0;
        if (csmin[i][0] > ctarget[i] || csmax[i][0] < ctarget[i]) empty = 1;   /* budget unreachable at the root */
    }
    int r = empty ? 0 : solve_rec(0, R);
    printf("%s %lld", r == 1 ? "FOUND" : (r == 2 ? "LIMIT" : "EXCLUDED"), nodes);
    for (int d = 0; d <= NC; d++) printf(" %lld", depthc[d]);
    printf(" | cons %d pruned %lld dyn %lld\n", ncons, cons_pruned, dyn_pruned);
    if (r == 1)
        for (int x = 0; x < NC; x++) {
            printf("%d %d", CC[x], CD[x]);
            for (int i = 0; i < 34; i++) printf(" %d", dom[x][34 * chosen[x] + i]);
            printf("\n");
        }
    return 0;
}
