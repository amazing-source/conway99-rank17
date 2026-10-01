/* Type-level test for the lattices L containing A7 with all 28 point pairs realised by roots (structures, 01/10).
   Each frame vector of cell {c,d} projects on span(A7) to t in A7^*; with tau = 4t in Z^8 (coordinates = the 8 points,
   point 7 = the omitted one), the options are  +-(3,3,-1^6) [c,d at 3],  +-(4 e_c - 4 e_d),  +-(2,2,6@z,-2^5).
   Necessary conditions used:  sum_U tau tau^T = 224 I - 28 J  (tight frame on the A7 span),
   star sums sum_{U in star(x)} tau_U = 0,  and for pairs the Cauchy-Schwarz window on the complementary part.
   Counts the type assignments of the 42 orbits (2 per cell) passing everything, up to a cap.                          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static int cc[21], cd[21];
static int opt[21][16][8], on2[21][16];       /* tau options per cell, squared norm (in 1/16 units) */
static int ch[42]; static long long nodes = 0, sols = 0, cap = 2000000000LL;
static int S[8][8];                            /* running sum of tau tau^T */
static int starsum[7][8], starleft[7];
static int maxrem[43][8];                      /* suffix max of tau_x^2 over remaining orbits */
static int ov(int x, int y) { return (cc[x]==cc[y]) + (cc[x]==cd[y]) + (cd[x]==cc[y]) + (cd[x]==cd[y]); }
static int pair_ok(int x, int a, int y, int b) {
    /* G = 1 - s - 2B for some B in {0,1,2};  16 G - tau.tau' = 16 (k,k'),  |(k,k')| <= |k||k'|, |k|^2 = 5 - |tau|^2/16 */
    int s = ov(x, y), d = 0;
    for (int i = 0; i < 8; i++) d += opt[x][a][i] * opt[y][b][i];
    double ka = 80 - on2[x][a], kb = 80 - on2[y][b];      /* 16 |k|^2 */
    double lim = sqrt(ka * kb) + 1e-9;
    for (int B = 0; B <= 2; B++) { double v = 16.0 * (1 - s - 2 * B) - d; if (fabs(v) <= lim) return 1; }
    return 0;
}
static void rec(int k) {
    nodes++;
    if (nodes > cap) return;
    if (k == 42) {
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (S[i][j] != (i == j ? 196 : -28)) return;
        sols++;
        if (sols <= 3) { printf("solution:"); for (int u = 0; u < 42; u++) printf(" %d", ch[u]); printf("\n"); }
        return;
    }
    int x = k / 2;
    for (int a = (k & 1) ? ch[k - 1] : 0; a < 16; a++) {      /* sisters unordered: second >= first */
        const int *t = opt[x][a]; int ok = 1;
        for (int i = 0; i < 8 && ok; i++) {
            int v = S[i][i] + t[i] * t[i];
            if (v > 196 || v + maxrem[k + 1][i] < 196) ok = 0;
        }
        if (!ok) continue;
        for (int j = 0; j < k && ok; j++) if (!pair_ok(x, a, j / 2, ch[j])) ok = 0;
        if (!ok) continue;
        /* off-diagonal Cauchy-Schwarz on the residual: |R_ij| <= sqrt(R_ii R_jj) at the end is exact; partial: skip */
        int c = cc[x], d = cd[x];
        for (int z = 0; z < 2 && ok; z++) {
            int l = z ? d : c;
            if (starleft[l] == 1) for (int i = 0; i < 8; i++) if (starsum[l][i] + t[i] != 0) { ok = 0; break; }
        }
        if (!ok) continue;
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) S[i][j] += t[i] * t[j];
        for (int i = 0; i < 8; i++) { starsum[c][i] += t[i]; starsum[d][i] += t[i]; }
        starleft[c]--; starleft[d]--; ch[k] = a;
        rec(k + 1);
        starleft[c]++; starleft[d]++;
        for (int i = 0; i < 8; i++) { starsum[c][i] -= t[i]; starsum[d][i] -= t[i]; }
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) S[i][j] -= t[i] * t[j];
    }
}
int main(void) {
    int x = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { cc[x] = c; cd[x] = d; x++; }
    for (x = 0; x < 21; x++) {
        int c = cc[x], d = cd[x], n = 0;
        for (int sg = 1; sg >= -1; sg -= 2) {
            for (int i = 0; i < 8; i++) opt[x][n][i] = sg * ((i == c || i == d) ? 3 : -1);
            n++;
            for (int i = 0; i < 8; i++) opt[x][n][i] = sg * (i == c ? 4 : (i == d ? -4 : 0));
            n++;
            for (int z = 0; z < 8; z++) {
                if (z == c || z == d) continue;
                for (int i = 0; i < 8; i++) opt[x][n][i] = sg * ((i == c || i == d) ? 2 : (i == z ? 6 : -2));
                n++;
            }
        }
        for (int a = 0; a < 16; a++) { int s = 0; for (int i = 0; i < 8; i++) s += opt[x][a][i] * opt[x][a][i]; on2[x][a] = s; }
    }
    for (int i = 0; i < 8; i++) maxrem[42][i] = 0;
    for (int k = 41; k >= 0; k--) for (int i = 0; i < 8; i++) {
        int m = 0; for (int a = 0; a < 16; a++) { int v = opt[k / 2][a][i] * opt[k / 2][a][i]; if (v > m) m = v; }
        maxrem[k][i] = maxrem[k + 1][i] + m; }
    for (int l = 0; l < 7; l++) starleft[l] = 12;
    rec(0);
    printf("nodes %lld type assignments %lld%s\n", nodes, sols, nodes > cap ? " (CAP)" : "");
    return 0;
}
