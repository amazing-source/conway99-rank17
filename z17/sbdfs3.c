/* Stage-B exact DFS with precomputed pairwise compatibility bitsets (structures, 28/09).
   Same rules and same input format as sbdfs.c (see there); used for heavy skeletons (after sbdfs.exe hits a small cap).
   Forward checking = AND of bitsets.  Star look-ahead for 1 and 2 open cells at a vertex, PSD residual, star closure.        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define NC 21
typedef uint64_t u64;
static int CC[NC], CD[NC], ndom[NC], nw[NC];
static signed char *dom[NC];
static double *proj[NC];
static u64 *kp[NC], *km[NC];
static u64 *compat[NC][NC];           /* compat[X][Y] : ndom[X] rows of nw[Y] words, X<Y */
static u64 *lev[NC + 1][NC];
static double Bm[17 * 15];
static long long nodes = 0, cap, depthc[NC + 1];
static int chosen[NC], rrt[NC][NC][2][2], lastk[7];
static u64 *htab; static int hcap = 0;

static unsigned long long ckey(const signed char *u, int plus, const int *v) {
    unsigned long long k = 0;
    for (int i = 16; i >= 0; i--) {
        int c = u ? (plus ? u[i] + u[17 + i] : u[i] - u[17 + i]) : v[i];
        k = k * 9ULL + (unsigned long long)(c + 4);
    }
    return k;
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
        u64 *src = lev[k][Y], *dst = lev[k + 1][Y];
        u64 acc = 0;
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
                    int q = 64 * w + bit;
                    const signed char *u = dom[Y] + 34 * q;
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
                int nB = 0;
                for (int w = 0; w < nw[Bc]; w++) nB += __builtin_popcountll(bB[w]);
                int hs = 1;
                while (hs < 2 * nB + 2) hs <<= 1;
                if (hs > hcap) { free(htab); htab = malloc(sizeof(u64) * (size_t)hs); hcap = hs; }
                memset(htab, 0, sizeof(u64) * (size_t)hs);
                u64 *kb = (CC[Bc] == v) ? kp[Bc] : km[Bc];
                for (int w = 0; w < nw[Bc]; w++) {
                    u64 x = bB[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        u64 key = kb[64 * w + bit] + 1ULL;
                        int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                        while (htab[pos] && htab[pos] != key) pos = (pos + 1) & (hs - 1);
                        htab[pos] = key;
                    }
                }
                for (int w = 0; w < nw[A]; w++) {
                    u64 x = bA[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        const signed char *u = dom[A] + 34 * (64 * w + bit);
                        int need[17], okr = 1;
                        for (int i = 0; i < 17; i++) {
                            int c = (CC[A] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                            need[i] = -s[i] - c;
                            if (need[i] < -4 || need[i] > 4) { okr = 0; break; }
                        }
                        int found = 0;
                        if (okr) {
                            u64 key = ckey(NULL, 0, need) + 1ULL;
                            int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                            while (htab[pos]) { if (htab[pos] == key) { found = 1; break; } pos = (pos + 1) & (hs - 1); }
                        }
                        if (!found) bA[w] &= ~(1ULL << bit);
                    }
                }
                if (!anybit(bA, nw[A])) return 0;
            }
        }
    }
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
    if (fread(Bm, sizeof(double), 17 * 15, f) != 17 * 15) return 1;
    int k = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { CC[k] = c; CD[k] = d; k++; }
    for (int x = 0; x < NC; x++) {
        if (fread(&ndom[x], 4, 1, f) != 1) return 1;
        int n = ndom[x] ? ndom[x] : 1;
        nw[x] = (n + 63) / 64;
        dom[x] = malloc(34 * (size_t)n);
        if (ndom[x] && fread(dom[x], 1, 34 * (size_t)ndom[x], f) != 34 * (size_t)ndom[x]) return 1;
        proj[x] = malloc(sizeof(double) * 30 * (size_t)n);
        kp[x] = malloc(sizeof(u64) * (size_t)n); km[x] = malloc(sizeof(u64) * (size_t)n);
        for (int t = 0; t < ndom[x]; t++) {
            kp[x][t] = ckey(dom[x] + 34 * t, 1, NULL); km[x][t] = ckey(dom[x] + 34 * t, 0, NULL);
            for (int s2 = 0; s2 < 2; s2++)
                for (int j = 0; j < 15; j++) {
                    double a = 0;
                    for (int i = 0; i < 17; i++) a += dom[x][34 * t + 17 * s2 + i] * Bm[i * 15 + j];
                    proj[x][30 * t + 15 * s2 + j] = a;
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
    if (!empty) {
        for (int X = 0; X < NC; X++)
            for (int Y = X + 1; Y < NC; Y++) {
                compat[X][Y] = calloc((size_t)ndom[X] * nw[Y], sizeof(u64));
                for (int t = 0; t < ndom[X]; t++) {
                    const signed char *a = dom[X] + 34 * t;
                    u64 *row = compat[X][Y] + (size_t)t * nw[Y];
                    for (int q = 0; q < ndom[Y]; q++) {
                        const signed char *b = dom[Y] + 34 * q;
                        int ok = 1;
                        for (int s = 0; s < 2 && ok; s++)
                            for (int tt = 0; tt < 2 && ok; tt++) {
                                int g = rrt[X][Y][s][tt];
                                for (int i = 0; i < 17; i++) g += a[17 * s + i] * b[17 * tt + i];
                                if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                            }
                        if (ok) row[q >> 6] |= 1ULL << (q & 63);
                    }
                }
            }
    }
    for (int v = 0; v < 7; v++) lastk[v] = -1;
    for (int x = 0; x < NC; x++) { lastk[CC[x]] = x; lastk[CD[x]] = x; }
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    int r = empty ? 0 : solve_rec(0, R);
    printf("%s %lld", r == 1 ? "FOUND" : (r == 2 ? "LIMIT" : "EXCLUDED"), nodes);
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
