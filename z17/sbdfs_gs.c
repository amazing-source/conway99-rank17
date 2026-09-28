/* Stage-B exact DFS for a general rank-17 unimodular V (structures, 28/09).  Rules as in SPEC_ENGINE.md / sbdfs3.c:
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
/* ---- Sol's graph-origin star rules (28/09, PROVEN from lambda = 1, mu = 2), option "star" ----
   Labels u = 2X + s (s = 0: r = e_c + e_d, s = 1: r = e_c - e_d). Local orientation of star z: sig_z(u) = (r_u)_z.
   (R1) within star z (both labels contain z): N_loc = sig_z(u) sig_z(v) (u,v); each u has #(+1) = #(-1) <= 1 there.
   (R2) u outside star z: N_loc = sig_z(v) (u,v) over v in star z; #(+1) = #(-1) <= 2.
   (R3) each u has exactly 10 labels v with (u,v) != 0.
   Upper bounds are pruned on the fly; equalities are checked once star z is complete; everything is checked at the leaf. */
static int use_star = 0;
static int lab_in[42][7], lab_sig[42][7];
static int cntOwn[42][7][2], cntOut[42][7][2], degc[42];
static int pstk[NC][200][3], pstn[NC];
static long long star_prunes = 0;
static int pairval(int X, int s, int tX, int Y, int t, int tY) {
    const signed char *a = dom[X] + 34 * tX + 17 * s;
    const int *gb = gdom[Y] + 34 * tY + 17 * t;
    int g = rrt[X][Y][s][t];
    for (int i = 0; i < 17; i++) g += a[i] * gb[i];
    return g / 4;                       /* g in {0,+-4} by the pairwise rule */
}
static void apply_pair(int u, int v, int p, int d) {
    degc[u] += d; degc[v] += d;
    for (int z = 0; z < 7; z++) {
        int iu = lab_in[u][z], iv = lab_in[v][z];
        if (iu && iv) { int loc = lab_sig[u][z] * lab_sig[v][z] * p, ix = loc > 0 ? 0 : 1; cntOwn[u][z][ix] += d; cntOwn[v][z][ix] += d; }
        else if (iu) { int loc = lab_sig[u][z] * p; cntOut[v][z][loc > 0 ? 0 : 1] += d; }
        else if (iv) { int loc = lab_sig[v][z] * p; cntOut[u][z][loc > 0 ? 0 : 1] += d; }
    }
}
static void star_remove(int k) {
    for (int i = pstn[k] - 1; i >= 0; i--) apply_pair(pstk[k][i][0], pstk[k][i][1], pstk[k][i][2], -1);
    pstn[k] = 0;
}
static int lab_ok_upper(int u) {
    if (degc[u] > 10) return 0;
    for (int z = 0; z < 7; z++) {
        if (lab_in[u][z]) { if (cntOwn[u][z][0] > 1 || cntOwn[u][z][1] > 1) return 0; }
        else { if (cntOut[u][z][0] > 2 || cntOut[u][z][1] > 2) return 0; }
    }
    return 1;
}
static int star_add(int k, int t) {          /* labels of cell k (entry t) against all assigned labels; 0 = violation */
    pstn[k] = 0;
    for (int s = 0; s < 2; s++) {
        int u = 2 * k + s;
        for (int Y = 0; Y < k; Y++)
            for (int q = 0; q < 2; q++) {
                int p = pairval(k, s, t, Y, q, chosen[Y]);
                if (p) { int v = 2 * Y + q; apply_pair(u, v, p, 1); pstk[k][pstn[k]][0] = u; pstk[k][pstn[k]][1] = v; pstk[k][pstn[k]][2] = p; pstn[k]++; }
            }
    }
    { int p = pairval(k, 0, t, k, 1, t); if (p) { apply_pair(2 * k, 2 * k + 1, p, 1); pstk[k][pstn[k]][0] = 2 * k; pstk[k][pstn[k]][1] = 2 * k + 1; pstk[k][pstn[k]][2] = p; pstn[k]++; } }
    int nass = 2 * (k + 1), rem = 2 * NC - nass;
    for (int u = 0; u < nass; u++) {
        if (!lab_ok_upper(u)) return 0;
        if (degc[u] + rem < 10) return 0;
    }
    for (int z = 0; z < 7; z++) {
        if (lastk[z] > k) continue;          /* star z complete: every assigned label's count into star z is final */
        for (int u = 0; u < nass; u++) {
            if (lab_in[u][z]) { if (cntOwn[u][z][0] != cntOwn[u][z][1]) return 0; }
            else { if (cntOut[u][z][0] != cntOut[u][z][1]) return 0; }
        }
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
    if (use_star) { if (!star_add(k, t)) { star_remove(k); star_prunes++; return 0; } }
    int rr = solve_rec(k + 1, R2);
    if (use_star && rr == 0) star_remove(k);
    return rr;
}

static int solve_rec(int k, double R[15][15]) {
    nodes++;
    depthc[k]++;
    if (nodes > cap) return 2;
    if (k == NC) {
        for (int i = 0; i < 15; i++)
            for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0;
        if (use_star)
            for (int u = 0; u < 2 * NC; u++) if (degc[u] != 10) return 0;
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
    for (int a = 2; a < argc; a++) if (!strcmp(argv[a], "star")) use_star = 1;
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
    for (int x = 0; x < NC; x++)
        for (int s2 = 0; s2 < 2; s2++) {
            int u = 2 * x + s2;
            lab_in[u][CC[x]] = 1; lab_sig[u][CC[x]] = 1;
            lab_in[u][CD[x]] = 1; lab_sig[u][CD[x]] = s2 ? -1 : 1;
        }
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    int r = empty ? 0 : solve_rec(0, R);
    printf("%s %lld", r == 1 ? "FOUND" : (r == 2 ? "LIMIT" : "EXCLUDED"), nodes);
    if (use_star) printf(" star_prunes %lld", star_prunes);
    for (int d = 0; d <= NC; d++) printf(" %lld", depthc[d]);
    printf("\n");
    if (r == 1)
        for (int x = 0; x < NC; x++) {
            printf("%d %d", CC[x], CD[x]);
            for (int i = 0; i < 34; i++) printf(" %d", dom[x][34 * chosen[x] + i]);
            printf("\n");
        }
    return 0;
}
