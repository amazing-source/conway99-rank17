/* evenframe.c -- exact search for the EVEN frame (Layer Q) inside a given lattice (structures, 01/10).
   Unknown: 42 vectors beta_U (U = (cell {c,d}, sister 0/1)), norm 5, in a lattice with Gram G (rank n <= 16), with
     (beta_U, beta_V) = 1 - s_UV - 2 B_UV,  s = |cell(U) cap cell(V)|,  B in {0,1,2}  (so the allowed values are
        s = 0: {1,-1,-3};  s = 1: {0,-2,-4};  s = 2 (sisters): {-1,-3,-5}),
     line counts  sum_{V in star(c), V != U} B_UV = 2 if c in cell(U), else 4      (this is B W = 4J - 2W),
     at most one B = 2 per row (the D-partner; exactly one follows from the counts and norm 5),
     star sums    sum_{V in star(c)} beta_V = 0,
     tight frame  sum_U beta_U beta_U^T = 14 Id.
   A solution is exactly a solution of Layer Q (B = X + 2D) whose frame lattice lies in the given lattice.
   Both sisters of a cell draw from the same candidate list (the constraint system is invariant under swapping the two
   sisters of any cell, so the second sister must have a larger index).  Global sign: the candidate lists are closed under
   negation and the first vector is taken with index in the first half of its list (list given as [X; -X]) when sym = 1.
   Input (binary LE): int32 n; int32 G[n*n]; int32 nodecap_millions; int32 sym; 21 cells (lexicographic (c<d), c,d in 0..6):
     int32 cnt; cnt * n int8.
   Output: EXCLUDED/FOUND/LIMIT nodes depth-profile; if FOUND the 42 vectors.                                           */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64;
#define NO 42
#define MAXN 16
static int n, G[MAXN][MAXN], sym;
static double Lt[MAXN][MAXN];
static int cnt[21], cc[21], cd[21];
static signed char *vec[21]; static int *gv[21]; static double *pr[21];
static int ord[NO], ocell[NO], osis[NO];              /* order of orbits; orbit -> cell, sister */
static int chosen[NO];
static u64 *live[NO + 1][NO]; static int nwd[21];
static int lineQ[NO][7], Dc[NO], remain[NO][7];      /* remain[U][c]: unchosen orbits of star(c) other than U */
static long long nodes = 0, cap, dprof[NO + 1];
static int starleft[7], maxdepth = 99;
static int incell(int x, int c) { return cc[x] == c || cd[x] == c; }
static int overlap(int x, int y) { return incell(x, cc[y]) + incell(x, cd[y]); }

static int psd_ok(double R[MAXN][MAXN]) {
    double L[MAXN][MAXN];
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= i; j++) {
            double s = R[i][j] + (i == j ? 1e-6 : 0.0);
            for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
            if (i == j) { if (s <= 0) return 0; L[i][i] = sqrt(s); }
            else L[i][j] = s / L[j][j];
        }
    return 1;
}
static inline int ipf(int x, int t, int y, int q) {
    const signed char *a = vec[x] + (size_t)t * n; const int *g = gv[y] + (size_t)q * n; int s = 0;
    for (int i = 0; i < n; i++) s += a[i] * g[i];
    return s;
}
static inline int qof(int s, int ip) {          /* B value or -1 */
    int b = 1 - s - ip;
    if (b == 0) return 0; if (b == 2) return 1; if (b == 4) return 2; return -1;
}
static int rec(int k, double R[MAXN][MAXN]);

static int try_vec(int k, int t, double R[MAXN][MAXN]) {
    int o = ord[k], x = ocell[o];
    int addQ[NO], nq = 0;
    /* counts against chosen orbits */
    int myline[7], myD = 0;
    for (int c = 0; c < 7; c++) myline[c] = 0;
    for (int j = 0; j < k; j++) {
        int V = ord[j], y = ocell[V];
        int q = qof(overlap(x, y), ipf(x, t, y, chosen[V]));
        if (q < 0) return 0;
        addQ[j] = q;
        if (q == 2) { if (Dc[V] || myD) return 0; myD = 1; }
        myline[cc[y]] += q; myline[cd[y]] += q;
    }
    for (int c = 0; c < 7; c++) {
        int tgt = incell(x, c) ? 2 : 4;
        int rem = remain[o][c];                       /* unchosen orbits of star(c) other than o (o itself not counted) */
        if (myline[c] > tgt || myline[c] + 2 * rem < tgt) return 0;
        if (rem && myD && myline[c] + rem < tgt && 0) return 0;
    }
    /* effect on chosen rows V: lines of cell(o) */
    for (int j = 0; j < k; j++) {
        int V = ord[j], y = ocell[V];
        for (int z = 0; z < 2; z++) {
            int c = z ? cd[x] : cc[x];
            int tgt = incell(y, c) ? 2 : 4;
            int cur = lineQ[V][c] + addQ[j];
            int rem = remain[V][c] - 1;
            if (cur > tgt || cur + 2 * rem < tgt) return 0;
        }
    }
    const double *p = pr[x] + (size_t)t * n;
    double R2[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j2 = 0; j2 < n; j2++) R2[i][j2] = R[i][j2] - p[i] * p[j2];
    for (int i = 0; i < n; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    /* commit */
    chosen[o] = t;
    for (int j = 0; j < k; j++) {
        int V = ord[j];
        lineQ[V][cc[x]] += addQ[j]; lineQ[V][cd[x]] += addQ[j];
        if (addQ[j] == 2) Dc[V] = 1;
    }
    for (int c = 0; c < 7; c++) lineQ[o][c] = myline[c];
    Dc[o] = myD;
    for (int U = 0; U < NO; U++) if (U != o) { remain[U][cc[x]]--; remain[U][cd[x]]--; }
    starleft[cc[x]]--; starleft[cd[x]]--;
    int ok = 1;
    /* forward checking */
    for (int j = k + 1; j < NO && ok; j++) {
        int V = ord[j], y = ocell[V], s = overlap(x, y);
        u64 *src = live[k][V], *dst = live[k + 1][V], acc = 0;
        for (int w = 0; w < nwd[y]; w++) {
            u64 xx = src[w], keep = 0;
            while (xx) {
                int bit = __builtin_ctzll(xx); xx &= xx - 1;
                int q = 64 * w + bit;
                int b = qof(s, ipf(x, t, y, q));
                if (b < 0) continue;
                if (b == 2 && myD && 0) continue;
                if (osis[V] == 1 && ocell[V] == x && osis[o] == 0 && q <= t) continue;   /* sister order */
                keep |= 1ULL << bit;
            }
            dst[w] = keep; acc |= keep;
        }
        if (!acc) ok = 0;
    }
    /* star sums: complete stars must vanish; stars with one orbit left force it */
    for (int z = 0; z < 2 && ok; z++) {
        int c = z ? cd[x] : cc[x];
        if (starleft[c] > 1) continue;
        int s[MAXN]; memset(s, 0, sizeof s);
        int last = -1;
        for (int U = 0; U < NO; U++) {
            if (!incell(ocell[U], c)) continue;
            int pos = -1;
            for (int j = 0; j <= k; j++) if (ord[j] == U) { pos = j; break; }
            if (pos < 0) { last = U; continue; }
            const signed char *v = vec[ocell[U]] + (size_t)chosen[U] * n;
            for (int i = 0; i < n; i++) s[i] += v[i];
        }
        if (starleft[c] == 0) { for (int i = 0; i < n; i++) if (s[i]) { ok = 0; break; } }
        else {
            int y = ocell[last]; u64 *b = live[k + 1][last]; int any = 0;
            for (int w = 0; w < nwd[y]; w++) {
                u64 xx = b[w];
                while (xx) {
                    int bit = __builtin_ctzll(xx); xx &= xx - 1;
                    const signed char *v = vec[y] + (size_t)(64 * w + bit) * n;
                    int eq = 1;
                    for (int i = 0; i < n; i++) if (v[i] != -s[i]) { eq = 0; break; }
                    if (!eq) b[w] &= ~(1ULL << bit); else any = 1;
                }
            }
            if (!any) ok = 0;
        }
    }
    int r = ok ? rec(k + 1, R2) : 0;
    if (r) return r;
    /* undo */
    for (int j = 0; j < k; j++) {
        int V = ord[j];
        lineQ[V][cc[x]] -= addQ[j]; lineQ[V][cd[x]] -= addQ[j];
        if (addQ[j] == 2) Dc[V] = 0;
    }
    Dc[o] = 0;
    for (int U = 0; U < NO; U++) if (U != o) { remain[U][cc[x]]++; remain[U][cd[x]]++; }
    starleft[cc[x]]++; starleft[cd[x]]++;
    return 0;
}

static int rec(int k, double R[MAXN][MAXN]) {
    nodes++; dprof[k]++;
    if (nodes > cap) return 2;
    if (k == maxdepth) return 0;
    if (k == NO) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (fabs(R[i][j]) > 1e-6) return 0;
        return 1;
    }
    int o = ord[k], x = ocell[o];
    u64 *b = live[k][o];
    for (int w = 0; w < nwd[x]; w++) {
        u64 xx = b[w];
        while (xx) {
            int bit = __builtin_ctzll(xx); xx &= xx - 1;
            int t = 64 * w + bit;
            if (k == 0 && sym && t >= cnt[x] / 2) continue;
            int r = try_vec(k, t, R);
            if (r) return r;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    if (!f) { puts("ERR open"); return 1; }
    int capm;
    if (fread(&n, 4, 1, f) != 1) return 1;
    int Gf[MAXN * MAXN];
    if (fread(Gf, 4, n * n, f) != (size_t)(n * n)) return 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) G[i][j] = Gf[i * n + j];
    if (fread(&capm, 4, 1, f) != 1 || fread(&sym, 4, 1, f) != 1) return 1;
    cap = (long long)capm * 1000000LL;
    { double L[MAXN][MAXN]; memset(L, 0, sizeof L);
      for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
          double s = G[i][j]; for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
          if (i == j) L[i][i] = sqrt(s); else L[i][j] = s / L[j][j]; }
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Lt[i][j] = L[j][i]; }
    int x = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { cc[x] = c; cd[x] = d; x++; }
    int empty = 0;
    for (x = 0; x < 21; x++) {
        if (fread(&cnt[x], 4, 1, f) != 1) return 1;
        vec[x] = malloc((size_t)cnt[x] * n + 1);
        if (cnt[x] && fread(vec[x], 1, (size_t)cnt[x] * n, f) != (size_t)cnt[x] * n) return 1;
        gv[x] = malloc(sizeof(int) * ((size_t)cnt[x] * n + 1)); pr[x] = malloc(sizeof(double) * ((size_t)cnt[x] * n + 1));
        for (int t = 0; t < cnt[x]; t++)
            for (int i = 0; i < n; i++) {
                int s = 0; double p = 0;
                for (int j = 0; j < n; j++) { s += G[i][j] * vec[x][(size_t)t * n + j]; p += Lt[i][j] * vec[x][(size_t)t * n + j]; }
                gv[x][(size_t)t * n + i] = s; pr[x][(size_t)t * n + i] = p;
            }
        nwd[x] = (cnt[x] + 63) / 64;
        if (cnt[x] < 2) empty = 1;
    }
    fclose(f);
    /* orbit order: stars in order: for each cell (lexicographic), sister 0 then sister 1 */
    for (int k = 0; k < NO; k++) { ord[k] = k; ocell[k] = k / 2; osis[k] = k & 1; }
    if (argc > 2) {            /* optional custom cell order */
        FILE *g = fopen(argv[2], "r"); int c;
        for (int k = 0; k < 21 && g && fscanf(g, "%d", &c) == 1; k++) { ord[2 * k] = 2 * c; ord[2 * k + 1] = 2 * c + 1; }
        if (g) fclose(g);
    }
    for (int L = 0; L <= NO; L++) for (int U = 0; U < NO; U++) live[L][U] = calloc(nwd[ocell[U]] + 1, sizeof(u64));
    for (int U = 0; U < NO; U++) for (int t = 0; t < cnt[ocell[U]]; t++) live[0][U][t >> 6] |= 1ULL << (t & 63);
    for (int U = 0; U < NO; U++) for (int c = 0; c < 7; c++) remain[U][c] = 12 - (incell(ocell[U], c) ? 1 : 0);
    for (int c = 0; c < 7; c++) starleft[c] = 12;
    double R[MAXN][MAXN];
    if (getenv("MAXDEPTH")) maxdepth = atoi(getenv("MAXDEPTH"));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) R[i][j] = (i == j) ? 14.0 : 0.0;
    int r = empty ? 0 : rec(0, R);
    printf("%s %lld |", r == 1 ? "FOUND" : (r == 2 ? "LIMIT" : "EXCLUDED"), nodes);
    for (int k = 0; k <= NO; k++) printf(" %lld", dprof[k]);
    printf("\n");
    if (r == 1)
        for (int U = 0; U < NO; U++) {
            printf("%d %d %d :", cc[ocell[U]], cd[ocell[U]], osis[U]);
            for (int i = 0; i < n; i++) printf(" %d", vec[ocell[U]][(size_t)chosen[U] * n + i]);
            printf("\n");
        }
    return 0;
}
