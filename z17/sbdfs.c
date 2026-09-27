/* Stage-B exact DFS for V = Z^17 (structures, 28/09).  Same rules as stageB.solve:
   cells in lexicographic order (c<d, 0-based), domains of ordered signed pairs (w+, w-) given by the caller (already filtered by
   (w+,w-) in {0,+-4}, M1, refined L1 and the sign normalisation);
   pairwise rule ((w_s^X, w_t^Y) + r_s^X . r_t^Y)/4 in {0,+-1} (forward checking);
   star closure  sum_{d>c}(w+ + w-) + sum_{d<c}(w+ - w-) = 0  when the last cell of a vertex is set, and look-ahead when a vertex
   has exactly one open cell (its a or b is then forced);
   frame: R = 28 I_15 - sum v v^T (v = projection of w on an orthonormal basis of W) must stay PSD: Cholesky of R + 1e-6 I.
   Input file (binary, little endian):
     int32 cap; double B[17*15] (orthonormal basis of W, row-major 17x15);
     for each of the 21 cells in order: int32 n; then n records of 34 int8 (w+[17], w-[17]).
   Output (stdout): "VERDICT nodes d0 d1 ... d21" and, if found, 21 lines "c d w+ w-".                                      */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NC 21
static int CC[NC], CD[NC];
static int ndom[NC];
static signed char *dom[NC];          /* n*34 */
static double *proj[NC];              /* n*30 : projections of w+ and w- (15 each) */
static unsigned long long *kp[NC], *km[NC]; /* keys of w+ + w- and w+ - w- */
static unsigned long long *htab; static int hcap = 0;
static double Bm[17 * 15];
static long long nodes = 0, cap;
static long long depthc[NC + 1];
static int chosen[NC];
static int rrt[NC][NC][2][2];
static int lastk[7];
/* per-level domain index lists */
static int *lev[NC + 1][NC];
static int levn[NC + 1][NC];

static int use2 = 1;

static int cmpu64(const void *a, const void *b) {
    unsigned long long x = *(const unsigned long long *)a, y = *(const unsigned long long *)b;
    return x < y ? -1 : (x > y);
}

/* base-9 key of a contribution vector (entries in [-4,4]): from a domain record (plus: w+ + w-, else w+ - w-) or from v */
static unsigned long long ckey(const signed char *u, int plus, const int *v) {
    unsigned long long k = 0;
    for (int i = 16; i >= 0; i--) {
        int c = u ? (plus ? u[i] + u[17 + i] : u[i] - u[17 + i]) : v[i];
        k = k * 9ULL + (unsigned long long)(c + 4);
    }
    return k;
}

static int dot17(const signed char *a, const signed char *b) {
    int s = 0;
    for (int i = 0; i < 17; i++) s += a[i] * b[i];
    return s;
}

static int psd_ok(double R[15][15]) {
    double L[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j <= i; j++) {
            double s = R[i][j] + (i == j ? 1e-6 : 0.0);
            for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
            if (i == j) {
                if (s <= 0) return 0;
                L[i][i] = sqrt(s);
            } else L[i][j] = s / L[j][j];
        }
    return 1;
}

static int star_ok(int v, int k) {
    /* all cells containing v are among order[0..k] */
    int s[17];
    memset(s, 0, sizeof s);
    for (int x = 0; x <= k; x++) {
        const signed char *w = dom[x] + 34 * chosen[x];
        if (CC[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] + w[17 + i];
        else if (CD[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] - w[17 + i];
    }
    for (int i = 0; i < 17; i++) if (s[i]) return 0;
    return 1;
}

static int solve_rec(int k, double R[15][15]);

static int try_entry(int k, int t, double R[15][15]) {
    int X = k;
    const signed char *w = dom[X] + 34 * t;
    const double *pv = proj[X] + 30 * t;
    double R2[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pv[i] * pv[j] - pv[15 + i] * pv[15 + j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    chosen[X] = t;
    for (int v = 0; v < 7; v++)
        if (lastk[v] == k && !star_ok(v, k)) return 0;
    /* forward checking into level k+1 */
    for (int Y = k + 1; Y < NC; Y++) {
        int m = 0;
        int *src = lev[k][Y], ns = levn[k][Y], *dst = lev[k + 1][Y];
        for (int q = 0; q < ns; q++) {
            const signed char *u = dom[Y] + 34 * src[q];
            int ok = 1;
            for (int s = 0; s < 2 && ok; s++)
                for (int tt = 0; tt < 2 && ok; tt++) {
                    int g = dot17(w + 17 * s, u + 17 * tt) + rrt[X][Y][s][tt];
                    if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                }
            if (ok) dst[m++] = src[q];
        }
        if (!m) return 0;
        levn[k + 1][Y] = m;
    }
    /* star look-ahead with two open cells: contrib(Y1) + contrib(Y2) = -s; keep only entries with a partner */
    if (use2) for (int v = 0; v < 7; v++) {
        int op[2], cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CC[Y] == v || CD[Y] == v) { if (cnt < 2) op[cnt] = Y; cnt++; }
        if (cnt != 2) continue;
        int s[17];
        memset(s, 0, sizeof s);
        for (int x = 0; x <= k; x++) {
            const signed char *ww = dom[x] + 34 * chosen[x];
            if (CC[x] == v) for (int i = 0; i < 17; i++) s[i] += ww[i] + ww[17 + i];
            else if (CD[x] == v) for (int i = 0; i < 17; i++) s[i] += ww[i] - ww[17 + i];
        }
        for (int side = 0; side < 2; side++) {
            int A = op[side], Bc = op[1 - side];
            int *arrB = lev[k + 1][Bc], nB = levn[k + 1][Bc];
            int hs = 1;
            while (hs < 2 * nB + 2) hs <<= 1;
            if (hs > hcap) { free(htab); htab = malloc(sizeof(unsigned long long) * (size_t)hs); hcap = hs; }
            memset(htab, 0, sizeof(unsigned long long) * (size_t)hs);
            unsigned long long *kb = (CC[Bc] == v) ? kp[Bc] : km[Bc];
            for (int q = 0; q < nB; q++) {
                unsigned long long key = kb[arrB[q]] + 1ULL;
                unsigned long long h = (key * 0x9E3779B97F4A7C15ULL) >> 20;
                int pos = (int)(h & (unsigned long long)(hs - 1));
                while (htab[pos] && htab[pos] != key) pos = (pos + 1) & (hs - 1);
                htab[pos] = key;
            }
            int *arrA = lev[k + 1][A], nA = levn[k + 1][A], m = 0;
            for (int q = 0; q < nA; q++) {
                int need[17], okr = 1;
                const signed char *u = dom[A] + 34 * arrA[q];
                for (int i = 0; i < 17; i++) {
                    int c = (CC[A] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                    need[i] = -s[i] - c;
                    if (need[i] < -4 || need[i] > 4) { okr = 0; break; }
                }
                if (!okr) continue;
                unsigned long long key = ckey(NULL, 0, need) + 1ULL;
                unsigned long long h = (key * 0x9E3779B97F4A7C15ULL) >> 20;
                int pos = (int)(h & (unsigned long long)(hs - 1)), found = 0;
                while (htab[pos]) { if (htab[pos] == key) { found = 1; break; } pos = (pos + 1) & (hs - 1); }
                if (found) arrA[m++] = arrA[q];
            }
            if (!m) return 0;
            levn[k + 1][A] = m;
        }
    }
    /* star look-ahead: vertex with exactly one open cell */
    for (int v = 0; v < 7; v++) {
        int open = -1, cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CC[Y] == v || CD[Y] == v) { open = Y; cnt++; }
        if (cnt != 1) continue;
        int s[17];
        memset(s, 0, sizeof s);
        for (int x = 0; x <= k; x++) {
            const signed char *ww = dom[x] + 34 * chosen[x];
            if (CC[x] == v) for (int i = 0; i < 17; i++) s[i] += ww[i] + ww[17 + i];
            else if (CD[x] == v) for (int i = 0; i < 17; i++) s[i] += ww[i] - ww[17 + i];
        }
        int *arr = lev[k + 1][open], n = levn[k + 1][open], m = 0;
        for (int q = 0; q < n; q++) {
            const signed char *u = dom[open] + 34 * arr[q];
            int ok = 1;
            for (int i = 0; i < 17 && ok; i++) {
                int c = (CC[open] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                if (c != -s[i]) ok = 0;
            }
            if (ok) arr[m++] = arr[q];
        }
        if (!m) return 0;
        levn[k + 1][open] = m;
    }
    return solve_rec(k + 1, R2);
}

static int solve_rec(int k, double R[15][15]) {
    nodes++;
    depthc[k]++;
    if (nodes > cap) return 2;
    if (k == NC) {
        for (int i = 0; i < 15; i++)
            for (int j = 0; j < 15; j++)
                if (fabs(R[i][j]) > 1e-6) return 0;
        return 1;
    }
    for (int q = 0; q < levn[k][k]; q++) {
        int r = try_entry(k, lev[k][k][q], R);
        if (r) return r;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc > 2) use2 = atoi(argv[2]);
    FILE *f = fopen(argv[1], "rb");
    if (!f) { printf("ERR open\n"); return 1; }
    int capi;
    if (fread(&capi, 4, 1, f) != 1) return 1;
    cap = capi;
    if (fread(Bm, sizeof(double), 17 * 15, f) != 17 * 15) return 1;
    int k = 0;
    for (int c = 0; c < 7; c++)
        for (int d = c + 1; d < 7; d++) { CC[k] = c; CD[k] = d; k++; }
    for (int x = 0; x < NC; x++) {
        if (fread(&ndom[x], 4, 1, f) != 1) return 1;
        dom[x] = malloc(34 * (size_t)(ndom[x] ? ndom[x] : 1));
        if (ndom[x] && fread(dom[x], 1, 34 * (size_t)ndom[x], f) != 34 * (size_t)ndom[x]) return 1;
        proj[x] = malloc(sizeof(double) * 30 * (size_t)(ndom[x] ? ndom[x] : 1));
        for (int t = 0; t < ndom[x]; t++)
            for (int s = 0; s < 2; s++)
                for (int j = 0; j < 15; j++) {
                    double a = 0;
                    for (int i = 0; i < 17; i++) a += dom[x][34 * t + 17 * s + i] * Bm[i * 15 + j];
                    proj[x][30 * t + 15 * s + j] = a;
                }
        kp[x] = malloc(sizeof(unsigned long long) * (size_t)(ndom[x] ? ndom[x] : 1));
        km[x] = malloc(sizeof(unsigned long long) * (size_t)(ndom[x] ? ndom[x] : 1));
        for (int t = 0; t < ndom[x]; t++) { kp[x][t] = ckey(dom[x] + 34 * t, 1, NULL); km[x][t] = ckey(dom[x] + 34 * t, 0, NULL); }
        for (int L = 0; L <= NC; L++) lev[L][x] = malloc(sizeof(int) * (size_t)(ndom[x] ? ndom[x] : 1));
        for (int t = 0; t < ndom[x]; t++) lev[0][x][t] = t;
        levn[0][x] = ndom[x];
    }
    fclose(f);
    for (int X = 0; X < NC; X++)
        for (int Y = 0; Y < NC; Y++)
            for (int s = 0; s < 2; s++)
                for (int t = 0; t < 2; t++) {
                    int r1[7] = {0}, r2[7] = {0};
                    r1[CC[X]] = 1; r1[CD[X]] = s ? -1 : 1;
                    r2[CC[Y]] = 1; r2[CD[Y]] = t ? -1 : 1;
                    int v = 0;
                    for (int i = 0; i < 7; i++) v += r1[i] * r2[i];
                    rrt[X][Y][s][t] = v;
                }
    for (int v = 0; v < 7; v++) lastk[v] = -1;
    for (int x = 0; x < NC; x++) { lastk[CC[x]] = x; lastk[CD[x]] = x; }
    double R[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    int empty = 0;
    for (int x = 0; x < NC; x++) if (!ndom[x]) empty = 1;
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
