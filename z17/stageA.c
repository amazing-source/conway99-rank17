/* Stage A for V = Z^17 in C (structures, 28/09).  Same enumeration as run_pair.py (Python reference):
   rows pi_1 = 0, pi_2..pi_7 chosen by partition refinement over the y-coordinate classes (one representative per orbit of the
   within-class permutation group), S_7 breaking N_12 = min, N_13 <= ... <= N_17; conditions: wt(pi_c + pi_d) in {2,6,10};
   even intersection with the odd supports of y1, y2; every class met is USABLE (norm-10 vectors w != +-w' of K with that odd
   support, not (3,1), with (w,w') in {0,+-4}); PE items (B, kmax) pruned during the tree (singletons also need k >= 1);
   last row covers every never-covered part.  Leaves are written as 17 columns (2 hex digits each), one per line.
   Input file (text): y1[17] y2[17] npe (B kmax single)*npe nclasses (size idx...)*nclasses
   Usage: stageA.exe INPUT OUTPUT                                                                                          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int Y1[17], Y2[17];
static int O1 = 0, O2 = 0;
static int npe, peB[4096], peK[4096], peS[4096];
static unsigned char *usab;   /* 0 unknown, 1 yes, 2 no, indexed by mask */
static int **amask; static int *namask;   /* per class: distinct abar masks */
static int use_m1 = 1;
static long long m1_prunes = 0;
static long long leaves = 0, usable_calls = 0;
static FILE *out;

static int pc(unsigned x) { return __builtin_popcount(x); }

static int usable(unsigned S) {
    if (usab[S]) return usab[S] == 1;
    usable_calls++;
    int idx[17], n = 0, rest[17], nr = 0;
    for (int i = 0; i < 17; i++) { if (S >> i & 1) idx[n++] = i; else rest[nr++] = i; }
    /* enumerate vectors */
    static signed char V[20000][17];
    int nv = 0;
    int ntw = (n == 10) ? 0 : (n == 6 ? 1 : (n == 2 ? 2 : -1));
    if (ntw < 0) { usab[S] = 2; return 0; }
    for (int sg = 0; sg < (1 << n); sg++) {
        signed char w[17];
        memset(w, 0, 17);
        for (int a = 0; a < n; a++) w[idx[a]] = (sg >> a & 1) ? -1 : 1;
        if (ntw == 0) {
            int d1 = 0, d2 = 0;
            for (int i = 0; i < 17; i++) { d1 += w[i] * Y1[i]; d2 += w[i] * Y2[i]; }
            if (!d1 && !d2) memcpy(V[nv++], w, 17);
        } else if (ntw == 1) {
            for (int a = 0; a < nr; a++)
                for (int s = -2; s <= 2; s += 4) {
                    w[rest[a]] = (signed char)s;
                    int d1 = 0, d2 = 0;
                    for (int i = 0; i < 17; i++) { d1 += w[i] * Y1[i]; d2 += w[i] * Y2[i]; }
                    if (!d1 && !d2) memcpy(V[nv++], w, 17);
                    w[rest[a]] = 0;
                }
        } else {
            for (int a = 0; a < nr; a++)
                for (int b = a + 1; b < nr; b++)
                    for (int s = -2; s <= 2; s += 4)
                        for (int t = -2; t <= 2; t += 4) {
                            w[rest[a]] = (signed char)s; w[rest[b]] = (signed char)t;
                            int d1 = 0, d2 = 0;
                            for (int i = 0; i < 17; i++) { d1 += w[i] * Y1[i]; d2 += w[i] * Y2[i]; }
                            if (!d1 && !d2) memcpy(V[nv++], w, 17);
                            w[rest[a]] = 0; w[rest[b]] = 0;
                        }
        }
    }
    int ok = 0;
    static unsigned char seen[1 << 17];
    static int list[1 << 17];
    int nl = 0;
    for (int p = 0; p < nv; p++)
        for (int q = 0; q < nv; q++) {
            if (p == q) continue;
            int g = 0;
            for (int i = 0; i < 17; i++) g += V[p][i] * V[q][i];
            if (g == 0 || g == 4 || g == -4) {
                ok = 1;
                int am = 0;
                for (int i = 0; i < 17; i++) { int a = (V[p][i] + V[q][i]) / 2; if (a & 1) am |= 1 << i; }
                if (!seen[am]) { seen[am] = 1; list[nl++] = am; }
            }
        }
    for (int t = 0; t < nl; t++) seen[list[t]] = 0;
    amask[S] = malloc(sizeof(int) * (size_t)(nl ? nl : 1));
    memcpy(amask[S], list, sizeof(int) * (size_t)nl);
    namask[S] = nl;
    usab[S] = ok ? 1 : 2;
    return ok;
}

/* partition: parts as arrays of coordinate indices */
typedef struct { int n; int idx[17]; } Part;

static unsigned rows[7];
static int nrows;

/* incremental M1: for every cell (c,d) with both rows known, some abar in its class must satisfy
   parity(abar & pi_e) == target_e(c) for all known rows e (target: block 0 -> e != 0; else e == c) */
static int m1_ok(void) {
    for (int c = 0; c < nrows; c++)
        for (int d = c + 1; d < nrows; d++) {
            unsigned S = rows[c] ^ rows[d];
            int *L = amask[S], n = namask[S], found = 0;
            for (int t = 0; t < n && !found; t++) {
                int good = 1;
                for (int e = 1; e < nrows && good; e++) {
                    int par = pc((unsigned)L[t] & rows[e]) & 1;
                    int tg = (c == 0) ? 1 : (e == c);
                    if (par != tg) good = 0;
                }
                if (good) found = 1;
            }
            if (!found) return 0;
        }
    return 1;
}

static int kfeasible(void) {
    int rem = 7 - nrows;
    for (int e = 0; e < npe; e++) {
        int w = 0;
        for (int r = 0; r < nrows; r++) if (pc(rows[r] & (unsigned)peB[e]) & 1) w++;
        int lo = peS[e] ? 1 : 0, ok = 0;
        for (int W = w; W <= w + rem && !ok; W++) {
            int k = W < 7 - W ? W : 7 - W;
            if (k >= lo && k <= peK[e]) ok = 1;
        }
        if (!ok) return 0;
    }
    return 1;
}

static void dfs(Part *parts, int np);

static int GOOD(int w) { return w == 2 || w == 6 || w == 10; }

/* choose next row: n_a from each part a (first n_a elements) */
static int selv[64];
static void choose(Part *parts, int np, int a, int wt, unsigned forced_mask) {
    if (wt > 10) return;
    if (a == np) {
        if (!GOOD(wt)) return;
        int c = nrows;                                /* index of new row */
        int w2 = c >= 2 ? pc(rows[1]) : 2;
        int wmin = (c <= 2) ? w2 : (w2 > pc(rows[c - 1]) ? w2 : pc(rows[c - 1]));
        if (c == 1) wmin = 2;
        if (wt < wmin) return;
        unsigned m = 0;
        for (int b = 0; b < np; b++) for (int j = 0; j < selv[b]; j++) m |= 1u << parts[b].idx[j];
        if ((pc(m & (unsigned)O1) & 1) || (pc(m & (unsigned)O2) & 1)) return;
        if (!usable(m)) return;
        for (int r = 1; r < c; r++) {
            int d = pc(m ^ rows[r]);
            if (!GOOD(d) || (c >= 2 && d < w2) || !usable(m ^ rows[r])) return;
        }
        rows[nrows++] = m;
        int m1 = !use_m1 || m1_ok();
        if (!m1) m1_prunes++;
        if (m1 && kfeasible()) {
            Part np2[64];
            int n2 = 0;
            for (int b = 0; b < np; b++) {
                if (selv[b]) { np2[n2].n = selv[b]; memcpy(np2[n2].idx, parts[b].idx, sizeof(int) * selv[b]); n2++; }
                if (selv[b] < parts[b].n) {
                    np2[n2].n = parts[b].n - selv[b];
                    memcpy(np2[n2].idx, parts[b].idx + selv[b], sizeof(int) * (parts[b].n - selv[b]));
                    n2++;
                }
            }
            dfs(np2, n2);
        }
        nrows--;
        return;
    }
    if (forced_mask >> a & 1) { selv[a] = parts[a].n; choose(parts, np, a + 1, wt + parts[a].n, forced_mask); selv[a] = 0; return; }
    for (int s = 0; s <= parts[a].n; s++) { selv[a] = s; choose(parts, np, a + 1, wt + s, forced_mask); }
    selv[a] = 0;
}

static void dfs(Part *parts, int np) {
    if (nrows == 7) {
        int cols[17];
        for (int i = 0; i < 17; i++) {
            int c = 0;
            for (int r = 0; r < 7; r++) if (rows[r] >> i & 1) c |= 1 << r;
            if (!c) return;
            cols[i] = c;
        }
        leaves++;
        for (int i = 0; i < 17; i++) fprintf(out, "%02x", cols[i]);
        fputc('\n', out);
        return;
    }
    unsigned forced = 0;
    if (nrows == 6) {
        unsigned cov = 0;
        for (int r = 0; r < 6; r++) cov |= rows[r];
        for (int a = 0; a < np; a++) if (!(cov >> parts[a].idx[0] & 1)) forced |= 1u << a;
    }
    /* selv is shared: save/restore across recursion levels */
    int save[64];
    memcpy(save, selv, sizeof save);
    memset(selv, 0, sizeof selv);
    choose(parts, np, 0, 0, forced);
    memcpy(selv, save, sizeof save);
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    out = fopen(argv[2], "w");
    for (int i = 0; i < 17; i++) if (fscanf(f, "%d", &Y1[i]) != 1) return 1;
    for (int i = 0; i < 17; i++) if (fscanf(f, "%d", &Y2[i]) != 1) return 1;
    for (int i = 0; i < 17; i++) { if (Y1[i] & 1) O1 |= 1 << i; if (Y2[i] & 1) O2 |= 1 << i; }
    if (fscanf(f, "%d", &npe) != 1) return 1;
    for (int e = 0; e < npe; e++) if (fscanf(f, "%d %d %d", &peB[e], &peK[e], &peS[e]) != 3) return 1;
    int ncl;
    if (fscanf(f, "%d", &ncl) != 1) return 1;
    Part parts[64];
    for (int a = 0; a < ncl; a++) {
        if (fscanf(f, "%d", &parts[a].n) != 1) return 1;
        for (int j = 0; j < parts[a].n; j++) if (fscanf(f, "%d", &parts[a].idx[j]) != 1) return 1;
    }
    fclose(f);
    usab = calloc(1u << 17, 1);
    amask = calloc(1u << 17, sizeof(int *)); namask = calloc(1u << 17, sizeof(int));
    if (argc > 3) use_m1 = atoi(argv[3]);
    rows[0] = 0; nrows = 1;
    dfs(parts, ncl);
    fclose(out);
    printf("leaves %lld usable_classes_computed %lld\n", leaves, usable_calls);
    return 0;
}
