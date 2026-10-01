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
static int dynstart = 12;
/* explicit symmetry elements: v -> v M (coefficient rows) with a permutation of the 7 lines; used for the lexmin test
   on the static prefix (cells in order, slot 0 then slot 1, values = candidate indices) */
#define GMAXCAP 20000
static int ngr = 0, GMAX = 3000; static short (*GM)[MAXN][MAXN]; static unsigned char (*GL)[7];
static int *gcinv; static short **gmap; static unsigned char asg[NO]; static long long canpruned = 0;
static int cellidx[7][7];
static int canon_ok(void) {
    for (int g = 1; g < ngr; g++) {
        const int *ci = gcinv + (size_t)g * 21;
        for (int c2 = 0; c2 < 21; c2++) {
            int c = ci[c2];
            if (!asg[2 * c]) break;
            const short *mp = gmap[g] + 0;
            int a = gmap[(size_t)g * 21 + c][chosen[2 * c]];
            (void)mp;
            if (asg[2 * c + 1]) {
                int b = gmap[(size_t)g * 21 + c][chosen[2 * c + 1]];
                int m0 = a < b ? a : b, m1 = a < b ? b : a;
                if (!asg[2 * c2]) break;
                if (m0 < chosen[2 * c2]) return 0;
                if (m0 > chosen[2 * c2]) break;
                if (!asg[2 * c2 + 1]) break;
                if (m1 < chosen[2 * c2 + 1]) return 0;
                if (m1 > chosen[2 * c2 + 1]) break;
            } else {
                if (!asg[2 * c2]) break;
                if (a < chosen[2 * c2]) return 0;
                break;
            }
        }
    }
    return 1;
}
static int starleft[7], maxdepth = 99, dbgdepth = -1, dbgleft = 0;
static signed char *QT[21][21];        /* QT[x][y][t*cnt[y]+q] = B value (0,1,2) or -1 */
static int ngen = 0; static short *gperm[21]; static signed char *gvec; static int *ggv; static int gnorm[512];
static unsigned char gact[NO + 1][512];
static long long sympruned = 0, budpruned = 0, cnpruned = 0;
static short cnsave[NO + 1][NO];
static signed char Bm[NO][NO]; static short CN[NO][NO];   /* B among chosen orbits; weighted common neighbours */
static int nbud = 0, btarget[256]; static short *bud[256][21]; static int brun[NO + 1][256];
static unsigned char *ccnt[NO];      /* per orbit: (NO+1) levels x cnt x 8 bytes (7 line counts + D flag) */
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
#define QF(x, t, y, q) (QT[x][y][(size_t)(t) * cnt[y] + (q)])
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
        int q = QF(x, t, y, chosen[V]);
        if (q < 0) return 0;
        addQ[j] = q;
        if (q == 2) { if (Dc[V] || myD) return 0; myD = 1; }
        myline[cc[y]] += q; myline[cd[y]] += q;
    }
    for (int c = 0; c < 7; c++) {
        int tgt = incell(x, c) ? 2 : 4;
        int rem = remain[o][c];                       /* unchosen orbits of star(c) other than o (o itself not counted) */
        if (myline[c] > tgt || myline[c] + rem + (myD ? 0 : 1) < tgt) return 0;
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
            if (cur > tgt || cur + rem + ((Dc[V] || addQ[j] == 2) ? 0 : 1) < tgt) return 0;
        }
    }
    /* common neighbours: (B^2)_UV = 4 - s_UV - B_UV for U != V, so every partial sum is bounded by it */
    {
        int cno[NO];
        for (int j = 0; j < k; j++) {
            int V = ord[j], y = ocell[V], sum = 0;
            for (int i = 0; i < k; i++) if (i != j && addQ[i]) sum += addQ[i] * Bm[ord[i]][V];
            cno[j] = sum;
            if (sum > 4 - overlap(x, y) - addQ[j]) { cnpruned++; return 0; }
        }
        for (int i = 0; i < k; i++) {
            if (!addQ[i]) continue;
            int U = ord[i];
            for (int j = i + 1; j < k; j++) {
                if (!addQ[j]) continue;
                int V = ord[j];
                if (CN[U][V] + addQ[i] * addQ[j] > 4 - overlap(ocell[U], ocell[V]) - Bm[U][V]) { cnpruned++; return 0; }
            }
        }
        for (int j = 0; j < k; j++) cnsave[k][j] = (short)cno[j];
    }
    if (ngr > 1 && k < dynstart) {
        int sv = chosen[o]; chosen[o] = t; asg[o] = 1;
        int okc = canon_ok();
        asg[o] = 0; chosen[o] = sv;
        if (!okc) { canpruned++; return 0; }
    }
    const double *p = pr[x] + (size_t)t * n;
    double R2[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j2 = 0; j2 < n; j2++) R2[i][j2] = R[i][j2] - p[i] * p[j2];
    for (int i = 0; i < n; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    /* commit */
    chosen[o] = t; asg[o] = 1;
    for (int j = 0; j < k; j++) {
        int V = ord[j];
        Bm[o][V] = Bm[V][o] = (signed char)addQ[j];
        CN[o][V] = CN[V][o] = cnsave[k][j];
    }
    for (int i = 0; i < k; i++) if (addQ[i]) for (int j = i + 1; j < k; j++) if (addQ[j]) {
        int U = ord[i], V = ord[j]; CN[U][V] += (short)(addQ[i] * addQ[j]); CN[V][U] = CN[U][V]; }
    for (int g = 0; g < ngen; g++) {
        unsigned char a = gact[k][g];
        if (a) { int d = 0; const signed char *v = vec[x] + (size_t)t * n; const int *gg = ggv + (size_t)g * n;
                 for (int i = 0; i < n; i++) d += v[i] * gg[i];
                 if (d) a = 0; }
        gact[k + 1][g] = a;
    }
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
    /* forward checking with per-candidate line counts */
    int satn = 0, satU[2 * NO], satc[2 * NO], dnew = -1;
    for (int j = 0; j < k; j++) {
        int V = ord[j], y = ocell[V];
        if (addQ[j] == 2) dnew = V;
        if (!addQ[j]) continue;
        for (int z = 0; z < 2; z++) {             /* rows whose count on a line of cell x changed at this step */
            int c = z ? cd[x] : cc[x];
            satU[satn] = V; satc[satn] = c; satn++;
        }
    }
    for (int j = k + 1; j < NO && ok; j++) {
        int V = ord[j], y = ocell[V], s = overlap(x, y);
        u64 *src = live[k][V], *dst = live[k + 1][V], acc = 0;
        unsigned char *c0 = ccnt[V] + (size_t)k * cnt[y] * 8, *c1 = ccnt[V] + (size_t)(k + 1) * cnt[y] * 8;
        int la = cc[x], lb = cd[x];
        int ta = incell(y, la) ? 2 : 4, tb = incell(y, lb) ? 2 : 4;
        int ra = remain[V][la], rb = remain[V][lb];
        int oa = cc[y], ob = cd[y];
        int roomA = (incell(x, oa) ? 2 : 4) - lineQ[o][oa], roomB = (incell(x, ob) ? 2 : 4) - lineQ[o][ob];
        for (int w = 0; w < nwd[y]; w++) {
            u64 xx = src[w], keep = 0;
            while (xx) {
                int bit = __builtin_ctzll(xx); xx &= xx - 1;
                int q = 64 * w + bit;
                int b = QF(x, t, y, q);
                if (b < 0) continue;
                if (ocell[V] == x && ((osis[o] == 0 && q <= t) || (osis[o] == 1 && q >= t))) continue;   /* sister order */
                const unsigned char *p0 = c0 + (size_t)q * 8; unsigned char *p1 = c1 + (size_t)q * 8;
                if (b == 2 && (p0[7] || Dc[o])) continue;
                if (b > roomA || b > roomB) continue;
                int na = p0[la] + b, nb = p0[lb] + b;
                { int dq = (p0[7] || b == 2) ? 0 : 1; if (na > ta || nb > tb || na + ra + dq < ta || nb + rb + dq < tb) continue; }
                memcpy(p1, p0, 8); p1[la] = (unsigned char)na; p1[lb] = (unsigned char)nb; if (b == 2) p1[7] = 1;
                /* rows saturated by this step, on a line of cell y */
                int bad = 0;
                for (int z = 0; z < satn && !bad; z++) {
                    if (satc[z] != oa && satc[z] != ob) continue;
                    int U = satU[z];
                    if (lineQ[U][satc[z]] + QF(ocell[U], chosen[U], y, q) > (incell(ocell[U], satc[z]) ? 2 : 4)) bad = 1;
                }
                if (!bad && dnew >= 0 && QF(ocell[dnew], chosen[dnew], y, q) == 2) bad = 1;
                if (bad) continue;
                keep |= 1ULL << bit;
            }
            dst[w] = keep; acc |= keep;
        }
        if (!acc) ok = 0;
    }
    /* budgets: sum over the 42 orbits of q_b(beta) = target (q_b = squared norm of the projection on a subspace) */
    if (ok && nbud) {
        for (int b_ = 0; b_ < nbud && ok; b_++) {
            int run = brun[k][b_] + bud[b_][x][t];
            brun[k + 1][b_] = run;
            int lo = run, hi = run;
            for (int j = k + 1; j < NO; j++) {
                int V = ord[j], y = ocell[V];
                const u64 *lv = live[k + 1][V]; const short *bc = bud[b_][y];
                int mn = 1 << 20, mx = -(1 << 20);
                for (int w = 0; w < nwd[y]; w++) {
                    u64 xx = lv[w];
                    while (xx) { int bit = __builtin_ctzll(xx); xx &= xx - 1; int v = bc[64 * w + bit]; if (v < mn) mn = v; if (v > mx) mx = v; }
                }
                lo += mn; hi += mx;
            }
            if (lo > btarget[b_] || hi < btarget[b_]) { ok = 0; budpruned++; }
        }
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
    asg[o] = 0;
    for (int i = 0; i < k; i++) if (addQ[i]) for (int j = i + 1; j < k; j++) if (addQ[j]) {
        int U = ord[i], V = ord[j]; CN[U][V] -= (short)(addQ[i] * addQ[j]); CN[V][U] = CN[U][V]; }
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
    int swapj = -1;
    if (k >= dynstart) {                 /* fail-first: smallest live domain next (valid: only static constraints from here) */
        int best = 1 << 30;
        for (int j = k; j < NO; j++) {
            int c = 0, V = ord[j];
            for (int w = 0; w < nwd[ocell[V]]; w++) c += __builtin_popcountll(live[k][V][w]);
            if (c < best) { best = c; swapj = j; }
        }
        int tmp = ord[k]; ord[k] = ord[swapj]; ord[swapj] = tmp;
    }
    int o = ord[k], x = ocell[o];
    u64 *b = live[k][o];
    if (dbgdepth == k && dbgleft > 0) {
        dbgleft--;
        fprintf(stderr, "live@%d:", k);
        for (int j = k; j < NO; j++) { int c = 0; for (int w = 0; w < nwd[ocell[ord[j]]]; w++) c += __builtin_popcountll(live[k][ord[j]][w]); fprintf(stderr, " %d", c); }
        fputc(10, stderr);
    }
    static short uf[NO + 1][2048];
    short *U_ = uf[k]; int anyg = 0;
    if (k < dynstart) for (int g = 0; g < ngen; g++) if (gact[k][g]) { anyg = 1; break; }
    if (anyg) {
        for (int q = 0; q < cnt[x]; q++) U_[q] = (short)q;
        for (int g = 0; g < ngen; g++) {
            if (!gact[k][g]) continue;
            const short *pm = gperm[x] + (size_t)g * cnt[x];
            for (int q = 0; q < cnt[x]; q++) {
                int a = q, c = pm[q];
                while (U_[a] != a) a = U_[a];
                while (U_[c] != c) c = U_[c];
                if (a < c) U_[c] = (short)a; else if (c < a) U_[a] = (short)c;
            }
        }
    }
    for (int w = 0; w < nwd[x]; w++) {
        u64 xx = b[w];
        while (xx) {
            int bit = __builtin_ctzll(xx); xx &= xx - 1;
            int t = 64 * w + bit;
            if (k == 0 && sym && t >= cnt[x] / 2) continue;
            if (anyg) { int a = t; while (U_[a] != a) a = U_[a]; if (a != t) { sympruned++; continue; } }
            int r = try_vec(k, t, R);
            if (r) return r;
        }
    }
    if (swapj >= 0) { int tmp = ord[k]; ord[k] = ord[swapj]; ord[swapj] = tmp; }
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
    if (fread(&ngen, 4, 1, f) != 1) ngen = 0;
    if (ngen > 512) { puts("ERR ngen"); return 1; }
    gvec = malloc((size_t)(ngen + 1) * n); ggv = malloc(sizeof(int) * (size_t)(ngen + 1) * n);
    if (ngen && fread(gvec, 1, (size_t)ngen * n, f) != (size_t)ngen * n) return 1;
    if (fread(&nbud, 4, 1, f) != 1) nbud = 0;
    if (nbud > 256) { puts("ERR nbud"); return 1; }
    for (int b_ = 0; b_ < nbud; b_++) {
        if (fread(&btarget[b_], 4, 1, f) != 1) return 1;
        for (int y = 0; y < 21; y++) { bud[b_][y] = malloc(2 * (size_t)cnt[y] + 2); if (fread(bud[b_][y], 2, cnt[y], f) != (size_t)cnt[y]) return 1; }
    }
    int ngg = 0; short (*GG)[MAXN][MAXN] = NULL; unsigned char (*GGL)[7] = NULL;
    if (fread(&ngg, 4, 1, f) != 1) ngg = 0;
    if (ngg) {
        GG = malloc(sizeof(*GG) * ngg); GGL = malloc(sizeof(*GGL) * ngg);
        for (int g = 0; g < ngg; g++) {
            signed char buf[MAXN * MAXN];
            if (fread(buf, 1, n * n, f) != (size_t)(n * n) || fread(GGL[g], 1, 7, f) != 7) return 1;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) GG[g][i][j] = buf[i * n + j];
        }
    }
    fclose(f);
    if (getenv("GMAX")) GMAX = atoi(getenv("GMAX"));
    if (GMAX > GMAXCAP) GMAX = GMAXCAP;
    { int z = 0; for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { cellidx[c][d] = cellidx[d][c] = z; z++; } }
    if (ngg) {
        GM = malloc(sizeof(*GM) * (GMAX + 1)); GL = malloc(sizeof(*GL) * (GMAX + 1));
        memset(GM[0], 0, sizeof(GM[0]));
        for (int i = 0; i < n; i++) GM[0][i][i] = 1;
        for (int c = 0; c < 7; c++) GL[0][c] = (unsigned char)c;
        ngr = 1;
        for (int h = 0; h < ngr && ngr < GMAX; h++)
            for (int g = 0; g < ngg && ngr < GMAX; g++) {
                short P[MAXN][MAXN]; unsigned char lp[7]; int big = 0;
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
                    int v = 0; for (int q = 0; q < n; q++) v += GM[h][i][q] * GG[g][q][j];
                    if (v > 30000 || v < -30000) big = 1;
                    P[i][j] = (short)v; }
                if (big) continue;
                for (int c = 0; c < 7; c++) lp[c] = GGL[g][GL[h][c]];
                int dup = 0;
                for (int e = 0; e < ngr && !dup; e++) if (!memcmp(GM[e], P, sizeof(short) * MAXN * MAXN) && !memcmp(GL[e], lp, 7)) dup = 1;
                if (dup) continue;
                memcpy(GM[ngr], P, sizeof P); memcpy(GL[ngr], lp, 7); ngr++;
            }
        /* candidate maps */
        gcinv = malloc(sizeof(int) * (size_t)ngr * 21); gmap = malloc(sizeof(short *) * (size_t)ngr * 21);
        int good = 1;
        for (int g = 0; g < ngr; g++) {
            int okg = 1;
            for (int c = 0; c < 21 && okg; c++) {
                int c2 = cellidx[GL[g][cc[c]]][GL[g][cd[c]]];
                gcinv[(size_t)g * 21 + c2] = c;
                short *mp = malloc(sizeof(short) * (cnt[c] + 1)); gmap[(size_t)g * 21 + c] = mp;
                for (int t = 0; t < cnt[c] && okg; t++) {
                    int img[MAXN]; const signed char *v = vec[c] + (size_t)t * n;
                    for (int j = 0; j < n; j++) { int a = 0; for (int i = 0; i < n; i++) a += v[i] * GM[g][i][j]; img[j] = a; }
                    int found = -1;
                    for (int q = 0; q < cnt[c2] && found < 0; q++) {
                        const signed char *u = vec[c2] + (size_t)q * n; int eq = 1;
                        for (int i = 0; i < n; i++) if (u[i] != img[i]) { eq = 0; break; }
                        if (eq) found = q;
                    }
                    if (found < 0) okg = 0; else mp[t] = (short)found;
                }
            }
            if (!okg) { good = 0; fprintf(stderr, "group element %d is not a symmetry of the instance%c", g, 10); break; }
        }
        if (!good) { puts("ERR group"); return 1; }
    }
    for (int g = 0; g < ngen; g++) {
        int nv = 0;
        for (int i = 0; i < n; i++) { int s2 = 0; for (int j = 0; j < n; j++) s2 += G[i][j] * gvec[(size_t)g * n + j]; ggv[(size_t)g * n + i] = s2; }
        for (int i = 0; i < n; i++) nv += gvec[(size_t)g * n + i] * ggv[(size_t)g * n + i];
        gnorm[g] = nv; gact[0][g] = 1;
    }
    for (x = 0; x < 21; x++) {
        if (cnt[x] > 2048) { puts("ERR cnt"); return 1; }
        gperm[x] = malloc(sizeof(short) * ((size_t)ngen * cnt[x] + 1));
        for (int g = 0; g < ngen; g++)
            for (int t = 0; t < cnt[x]; t++) {
                int d = 0, img[MAXN]; const signed char *v = vec[x] + (size_t)t * n;
                for (int i = 0; i < n; i++) d += v[i] * ggv[(size_t)g * n + i];
                int m = 2 * d / gnorm[g];
                for (int i = 0; i < n; i++) img[i] = v[i] - m * gvec[(size_t)g * n + i];
                int found = -1;
                for (int q = 0; q < cnt[x] && found < 0; q++) {
                    const signed char *u = vec[x] + (size_t)q * n; int eq = 1;
                    for (int i = 0; i < n; i++) if (u[i] != img[i]) { eq = 0; break; }
                    if (eq) found = q;
                }
                if (found < 0) { printf("ERR generator %d does not preserve cell %d%c", g, x, 10); return 1; }
                gperm[x][(size_t)g * cnt[x] + t] = (short)found;
            }
    }
    for (x = 0; x < 21; x++)
        for (int y = 0; y < 21; y++) {
            QT[x][y] = malloc((size_t)cnt[x] * cnt[y] + 1);
            int s2 = overlap(x, y);
            for (int t = 0; t < cnt[x]; t++)
                for (int q = 0; q < cnt[y]; q++) QT[x][y][(size_t)t * cnt[y] + q] = (signed char)qof(s2, ipf(x, t, y, q));
        }
    /* orbit order: stars in order: for each cell (lexicographic), sister 0 then sister 1 */
    for (int k = 0; k < NO; k++) { ord[k] = k; ocell[k] = k / 2; osis[k] = k & 1; }
    if (argc > 2) {            /* optional custom cell order */
        FILE *g = fopen(argv[2], "r"); int c;
        for (int k = 0; k < 21 && g && fscanf(g, "%d", &c) == 1; k++) { ord[2 * k] = 2 * c; ord[2 * k + 1] = 2 * c + 1; }
        if (g) fclose(g);
    }
    for (int L = 0; L <= NO; L++) for (int U = 0; U < NO; U++) live[L][U] = calloc(nwd[ocell[U]] + 1, sizeof(u64));
    for (int U = 0; U < NO; U++) ccnt[U] = calloc((size_t)(NO + 1) * cnt[ocell[U]] * 8 + 8, 1);
    for (int U = 0; U < NO; U++) for (int t = 0; t < cnt[ocell[U]]; t++) live[0][U][t >> 6] |= 1ULL << (t & 63);
    for (int U = 0; U < NO; U++) for (int c = 0; c < 7; c++) remain[U][c] = 12 - (incell(ocell[U], c) ? 1 : 0);
    for (int c = 0; c < 7; c++) starleft[c] = 12;
    double R[MAXN][MAXN];
    if (getenv("MAXDEPTH")) maxdepth = atoi(getenv("MAXDEPTH"));
    if (getenv("DYNSTART")) dynstart = atoi(getenv("DYNSTART"));
    if (getenv("DBGDEPTH")) { dbgdepth = atoi(getenv("DBGDEPTH")); dbgleft = 6; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) R[i][j] = (i == j) ? 14.0 : 0.0;
    for (int b_ = 0; b_ < nbud; b_++) {          /* static root check */
        int lo = 0, hi = 0;
        for (int U = 0; U < NO; U++) { int y = ocell[U], mn = 1 << 20, mx = -(1 << 20);
            for (int q = 0; q < cnt[y]; q++) { int v = bud[b_][y][q]; if (v < mn) mn = v; if (v > mx) mx = v; }
            lo += mn; hi += mx; }
        if (lo > btarget[b_] || hi < btarget[b_]) { empty = 1; budpruned++; }
        brun[0][b_] = 0;
    }
    int r = empty ? 0 : rec(0, R);
    printf("%s %lld sym %lld bud %lld cn %lld can %lld grp %d |", r == 1 ? "FOUND" : (r == 2 ? "LIMIT" : "EXCLUDED"), nodes, sympruned, budpruned, cnpruned, canpruned, ngr);
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
