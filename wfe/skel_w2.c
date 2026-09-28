/* skel_w2.c -- self-contained skeleton enumeration with witness budgets (structures, 28/09 night).
   Same search as skel_w.c (C1's skeleton rules + H-reduction + S7 base-line reduction + witness o-rules + unit rule), but
   the M1-aware budget tables are built here from the raw filtered shell:
     valid pairs of a class = unordered {w, w'} of shell vectors of the class with (w, w') in {0, +-4};
     a-masks m = ((w +- w')/2)^T GK mod 2 (two per pair); contribution for witness pair (t,t') =
     a_t(w)a_t'(w) + a_t(w')a_t'(w')  (the same for all 8 sign/order variants).
   Phase 1 (sampling): Knuth walks; at a walk's leaf level every M1-valid leaf is evaluated over ALL witness pairs (on the
   fly from the pair lists); the killing pairs of the sampled leaves are counted.  Phase 2: greedy double cover + most
   frequent -> P (<= PMAX pairs), int8 tables per (class, mask).  Phase 3: Knuth estimate or exhaustive run (sharded).
   Survivors (M1 + budgets over P) are printed for the exact check in Python (all pairs, then the vector level).
   Input (binary LE): int32 nvec, nw, nU, nun; int32 GK[225]; int32 target[nw(nw+1)/2] (pairs t<=t' row-major);
     int8 kall[nw]; per vector: int8 x[15], int8 val[nw], int32 cls_index (U order, sorted ascending);
     per class (U order): int32 cls, int32 Fm, int8 orbmin, uint64 wpar, uint32 upar.
   Usage: skel_w2.exe input mode(est|run) nshards shard survivors_out walks seed [pmax]                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
static int nvec, nw, nU, nun, NPW, nP, nwd, pmax = 40;
static int GK[15][15];
static int *targetAll, *target, *Plist;
static signed char *kall, *X, *VAL;
static int *vcls, *vstart, *vend;
static int *cls, *Fm; static signed char *orbmin; static u64 *wpar; static uint32_t *upar;
static int idxof[1 << 15];
/* pairs */
static int *pst, *pen, *PI, *PJ; static int16_t *PM1, *PM2; static long long npairs = 0;
/* per class unique masks and tables */
static int *mst, *men; static int16_t *UM; static int *R1, *R2;   /* row index of each pair's two masks */
static signed char *TMN, *TMX;                                   /* rows x nP */
static int *TPa, *TPb;                                           /* witness indices of all pairs */
static u64 *nbr; static u64 *candl[7];
static int FEAS[64][7][7];
static long long nodes = 0, leaves = 0, m1dead = 0, bdead = 0, surv = 0;
static int chosen[6];
static int nshards = 1, shard = 0;
static FILE *fout;
static int ocnt[7][64], ucnt[7][32];
static unsigned char spanmark[1 << 15];
static int spanel[64], nspan = 1;
static long long *killcnt; static int *lo, *hi, *cmn, *cmx;
static unsigned long long rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static int pc(u64 x) { return __builtin_popcountll(x); }
static int patof(int m, const int *G) { int p = 0; for (int e = 0; e < 6; e++) p |= (__builtin_popcount(m & G[e]) & 1) << e; return p; }

/* table-based leaf check over P */
static int leaf_check(void) {
    int G[6];
    for (int e = 0; e < 6; e++) G[e] = cls[chosen[e]];
    for (int p = 0; p < nP; p++) { lo[p] = 0; hi[p] = 0; }
    for (int c = -1; c < 6; c++)
        for (int d = c + 1; d < 6; d++) {
            int C = (c < 0) ? G[d] : (G[c] ^ G[d]);
            int i = idxof[C];
            if (i < 0) return 0;
            int want = (c < 0) ? 63 : (1 << c), any = 0;
            for (int p = 0; p < nP; p++) { cmn[p] = 1 << 20; cmx[p] = -(1 << 20); }
            for (int q = mst[i]; q < men[i]; q++) {
                if (patof(UM[q], G) != want) continue;
                any = 1;
                const signed char *a = TMN + (size_t)q * nP, *b = TMX + (size_t)q * nP;
                for (int p = 0; p < nP; p++) { if (a[p] < cmn[p]) cmn[p] = a[p]; if (b[p] > cmx[p]) cmx[p] = b[p]; }
            }
            if (!any) { m1dead++; return 0; }
            for (int p = 0; p < nP; p++) { lo[p] += cmn[p]; hi[p] += cmx[p]; }
        }
    for (int p = 0; p < nP; p++) if (lo[p] > target[p] || hi[p] < target[p]) { bdead++; return 0; }
    return 1;
}

/* M1 only (for sampling before tables exist): every cell has a matching mask */
static int leaf_m1(void) {
    int G[6];
    for (int e = 0; e < 6; e++) G[e] = cls[chosen[e]];
    for (int c = -1; c < 6; c++)
        for (int d = c + 1; d < 6; d++) {
            int C = (c < 0) ? G[d] : (G[c] ^ G[d]);
            int i = idxof[C], want = (c < 0) ? 63 : (1 << c), any = 0;
            for (int q = mst[i]; q < men[i] && !any; q++) if (patof(UM[q], G) == want) any = 1;
            if (!any) return 0;
        }
    return 1;
}

/* all-pair evaluation of the current leaf from the pair lists: counts killing witness pairs */
static int *alo, *ahi, *acmn, *acmx;
static int leaf_allpairs(void) {
    int G[6];
    for (int e = 0; e < 6; e++) G[e] = cls[chosen[e]];
    for (int p = 0; p < NPW; p++) { alo[p] = 0; ahi[p] = 0; }
    for (int c = -1; c < 6; c++)
        for (int d = c + 1; d < 6; d++) {
            int C = (c < 0) ? G[d] : (G[c] ^ G[d]);
            int i = idxof[C], want = (c < 0) ? 63 : (1 << c), any = 0;
            for (int p = 0; p < NPW; p++) { acmn[p] = 1 << 20; acmx[p] = -(1 << 20); }
            for (long long q = pst[i]; q < pen[i]; q++) {
                if (patof(PM1[q], G) != want && patof(PM2[q], G) != want) continue;
                any = 1;
                const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
                for (int p = 0; p < NPW; p++) {
                    int v = va[TPa[p]] * va[TPb[p]] + vb[TPa[p]] * vb[TPb[p]];
                    if (v < acmn[p]) acmn[p] = v;
                    if (v > acmx[p]) acmx[p] = v;
                }
            }
            if (!any) return -1;
            for (int p = 0; p < NPW; p++) { alo[p] += acmn[p]; ahi[p] += acmx[p]; }
        }
    int nk = 0;
    for (int p = 0; p < NPW; p++) if (alo[p] > targetAll[p] || ahi[p] < targetAll[p]) { killcnt[p]++; nk++; }
    return nk;
}

static void rec(int depth) {
    nodes++;
    u64 *cand = candl[depth];
    int need = 5 - depth;
    for (int w = 0; w < nwd; w++) {
        u64 x = cand[w];
        while (x) {
            int bit = __builtin_ctzll(x); x &= x - 1;
            int i = 64 * w + bit, g = cls[i];
            if (spanmark[g]) continue;
            if (depth == 0) {
                if (!orbmin[i]) continue;
                if (nshards > 1 && (i % nshards) != shard) continue;
            } else {
                int bad = 0;
                for (int j = 0; j < depth && !bad; j++) if (idxof[cls[chosen[j]] ^ g] < chosen[0]) bad = 1;
                if (bad) continue;
            }
            int r = 5 - depth, bad = 0;
            for (int t = 0; t < nw; t++) {
                int o = ocnt[depth][t] + (int)((wpar[i] >> t) & 1);
                ocnt[depth + 1][t] = o;
                if (!FEAS[t][o][r]) { bad = 1; break; }
            }
            if (bad) continue;
            for (int u = 0; u < nun; u++) ucnt[depth + 1][u] = ucnt[depth][u] + (int)((upar[i] >> u) & 1);
            chosen[depth] = i;
            if (depth == 5) {
                int uok = 1;
                for (int u = 0; u < nun; u++) if (ucnt[6][u] < 1) uok = 0;
                if (!uok) continue;
                leaves++;
                if (leaf_check()) {
                    surv++;
                    if (fout) { for (int e = 0; e < 6; e++) fprintf(fout, "%d ", cls[chosen[e]]); fputc(10, fout); fflush(fout); }
                }
                continue;
            }
            u64 *c2 = candl[depth + 1], *nb = nbr + (size_t)i * nwd;
            int cnt = 0;
            for (int k = 0; k < nwd; k++) { c2[k] = cand[k] & nb[k]; cnt += pc(c2[k]); }
            if (cnt < need) continue;
            int old = nspan;
            for (int s = 0; s < old; s++) { int y = spanel[s] ^ g; spanel[nspan++] = y; spanmark[y] = 1; }
            rec(depth + 1);
            for (int s = old; s < nspan; s++) spanmark[spanel[s]] = 0;
            nspan = old;
        }
    }
}

static int kids[1 << 15];
static double est_nodes, est_leaves, est_m1, est_surv;
static long long samp_leaves = 0, samp_killed = 0, samp_m1 = 0;

/* mode 0: sampling (all-pair evaluation of one random M1-valid leaf per walk); mode 1: estimate with tables */
static void walk(int mode) {
    double mult = 1.0;
    for (int s = 1; s < nspan; s++) spanmark[spanel[s]] = 0;
    nspan = 1; spanel[0] = 0; spanmark[0] = 1;
    for (int depth = 0; depth < 6; depth++) {
        est_nodes += mult;
        u64 *cand = candl[depth];
        int nk = 0;
        for (int w = 0; w < nwd; w++) {
            u64 x = cand[w];
            while (x) {
                int bit = __builtin_ctzll(x); x &= x - 1;
                int i = 64 * w + bit, g = cls[i];
                if (spanmark[g]) continue;
                if (depth == 0) { if (!orbmin[i]) continue; }
                else {
                    int bad = 0;
                    for (int j = 0; j < depth && !bad; j++) if (idxof[cls[chosen[j]] ^ g] < chosen[0]) bad = 1;
                    if (bad) continue;
                }
                int r = 5 - depth, bad = 0;
                for (int t = 0; t < nw; t++) if (!FEAS[t][ocnt[depth][t] + (int)((wpar[i] >> t) & 1)][r]) { bad = 1; break; }
                if (bad) continue;
                if (depth == 5) {
                    int uok = 1;
                    for (int u = 0; u < nun; u++) if (ucnt[depth][u] + (int)((upar[i] >> u) & 1) < 1) uok = 0;
                    if (!uok) continue;
                } else {
                    u64 *nb = nbr + (size_t)i * nwd; int cnt = 0;
                    for (int k = 0; k < nwd; k++) cnt += pc(cand[k] & nb[k]);
                    if (cnt < 5 - depth) continue;
                }
                kids[nk++] = i;
            }
        }
        if (!nk) return;
        if (depth == 5) {
            est_leaves += mult * nk;
            if (mode == 0) {
                int nm1 = 0, pick = -1;
                for (int q = 0; q < nk; q++) {
                    chosen[5] = kids[q];
                    if (leaf_m1()) { nm1++; if (pick < 0 || (xr() % nm1) == 0) pick = kids[q]; }
                }
                est_m1 += mult * nm1;
                if (pick >= 0) {
                    chosen[5] = pick;
                    samp_leaves++;
                    int r = leaf_allpairs();
                    if (r > 0) samp_killed++;
                    if (r < 0) samp_m1++;
                }
            } else {
                int ns = 0;
                for (int q = 0; q < nk; q++) { chosen[5] = kids[q]; if (leaf_check()) ns++; }
                est_surv += mult * ns;
            }
            return;
        }
        int i = kids[xr() % nk];
        mult *= nk;
        chosen[depth] = i;
        int g = cls[i];
        for (int t = 0; t < nw; t++) ocnt[depth + 1][t] = ocnt[depth][t] + (int)((wpar[i] >> t) & 1);
        for (int u = 0; u < nun; u++) ucnt[depth + 1][u] = ucnt[depth][u] + (int)((upar[i] >> u) & 1);
        u64 *c2 = candl[depth + 1], *nb = nbr + (size_t)i * nwd;
        for (int k = 0; k < nwd; k++) c2[k] = cand[k] & nb[k];
        int old = nspan;
        for (int s = 0; s < old; s++) { int y = spanel[s] ^ g; spanel[nspan++] = y; spanmark[y] = 1; }
    }
}

static int cmp16(const void *a, const void *b) { return (int)*(const int16_t *)a - (int)*(const int16_t *)b; }

int main(int argc, char **argv) {
    if (argc < 8) { puts("usage: skel_w2.exe input est|run nshards shard survivors walks seed [pmax]"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { puts("ERR open"); return 1; }
    int hdr[4];
    if (fread(hdr, 4, 4, f) != 4) return 1;
    nvec = hdr[0]; nw = hdr[1]; nU = hdr[2]; nun = hdr[3];
    if (nw > 64 || nun > 32) { puts("ERR sizes"); return 1; }
    if (fread(GK, 4, 225, f) != 225) return 1;
    NPW = nw * (nw + 1) / 2;
    targetAll = malloc(sizeof(int) * (NPW + 1)); TPa = malloc(sizeof(int) * (NPW + 1)); TPb = malloc(sizeof(int) * (NPW + 1));
    if (fread(targetAll, 4, NPW, f) != (size_t)NPW) return 1;
    { int p = 0; for (int t = 0; t < nw; t++) for (int t2 = t; t2 < nw; t2++) { TPa[p] = t; TPb[p] = t2; p++; } }
    kall = malloc(nw + 1);
    if (fread(kall, 1, nw, f) != (size_t)nw) return 1;
    X = malloc((size_t)nvec * 15); VAL = malloc((size_t)nvec * nw + 1); vcls = malloc(sizeof(int) * (size_t)nvec);
    for (int v = 0; v < nvec; v++)
        if (fread(X + (size_t)v * 15, 1, 15, f) != 15 || (nw && fread(VAL + (size_t)v * nw, 1, nw, f) != (size_t)nw) ||
            fread(vcls + v, 4, 1, f) != 1) { puts("ERR vec"); return 1; }
    cls = malloc(sizeof(int) * nU); Fm = malloc(sizeof(int) * nU); orbmin = malloc(nU); wpar = malloc(8 * (size_t)nU);
    upar = malloc(4 * (size_t)nU);
    for (int i = 0; i < (1 << 15); i++) idxof[i] = -1;
    for (int i = 0; i < nU; i++) {
        if (fread(&cls[i], 4, 1, f) != 1 || fread(&Fm[i], 4, 1, f) != 1 || fread(&orbmin[i], 1, 1, f) != 1 ||
            fread(&wpar[i], 8, 1, f) != 1 || fread(&upar[i], 4, 1, f) != 1) { puts("ERR class"); return 1; }
        idxof[cls[i]] = i;
    }
    fclose(f);
    nshards = atoi(argv[3]); shard = atoi(argv[4]);
    int walks = atoi(argv[6]);
    rng ^= (unsigned long long)atoll(argv[7]) * 0x9E3779B97F4A7C15ULL;
    if (argc > 8) pmax = atoi(argv[8]);
    /* vector ranges per class */
    vstart = malloc(sizeof(int) * nU); vend = malloc(sizeof(int) * nU);
    for (int i = 0; i < nU; i++) { vstart[i] = 0; vend[i] = 0; }
    for (int v = 0; v < nvec; v++) { int c = vcls[v]; if (v == 0 || vcls[v - 1] != c) vstart[c] = v; vend[c] = v + 1; }
    /* G w per vector */
    int *GW = malloc(sizeof(int) * 15 * (size_t)nvec);
    for (int v = 0; v < nvec; v++)
        for (int k = 0; k < 15; k++) { int s = 0; for (int j = 0; j < 15; j++) s += GK[k][j] * X[(size_t)v * 15 + j]; GW[(size_t)v * 15 + k] = s; }
    /* pairs per class */
    long long cap = 1 << 20;
    PI = malloc(sizeof(int) * cap); PJ = malloc(sizeof(int) * cap); PM1 = malloc(2 * cap); PM2 = malloc(2 * cap);
    pst = malloc(sizeof(int) * nU); pen = malloc(sizeof(int) * nU);
    for (int c = 0; c < nU; c++) {
        pst[c] = (int)npairs;
        for (int a = vstart[c]; a < vend[c]; a++)
            for (int b = a + 1; b < vend[c]; b++) {
                int ip = 0;
                for (int k = 0; k < 15; k++) ip += X[(size_t)a * 15 + k] * GW[(size_t)b * 15 + k];
                if (!(ip == 0 || ip == 4 || ip == -4)) continue;
                if (npairs == cap) {
                    cap *= 2;
                    PI = realloc(PI, sizeof(int) * cap); PJ = realloc(PJ, sizeof(int) * cap);
                    PM1 = realloc(PM1, 2 * cap); PM2 = realloc(PM2, 2 * cap);
                }
                int m1 = 0, m2 = 0;
                for (int k = 0; k < 15; k++) {
                    int s = GW[(size_t)a * 15 + k] + GW[(size_t)b * 15 + k], d = GW[(size_t)a * 15 + k] - GW[(size_t)b * 15 + k];
                    m1 |= ((s / 2) & 1) << k; m2 |= ((d / 2) & 1) << k;
                }
                PI[npairs] = a; PJ[npairs] = b; PM1[npairs] = (int16_t)m1; PM2[npairs] = (int16_t)m2; npairs++;
            }
        pen[c] = (int)npairs;
        if (npairs > 2000000000LL) { puts("ERR too many pairs"); return 1; }
    }
    /* unique masks per class */
    mst = malloc(sizeof(int) * nU); men = malloc(sizeof(int) * nU);
    UM = malloc(2 * (size_t)(2 * npairs + 1));
    int16_t *tmp = malloc(2 * (size_t)(2 * npairs + 1));
    long long nrows = 0;
    for (int c = 0; c < nU; c++) {
        int n = 0;
        for (long long q = pst[c]; q < pen[c]; q++) { tmp[n++] = PM1[q]; tmp[n++] = PM2[q]; }
        qsort(tmp, n, 2, cmp16);
        mst[c] = (int)nrows;
        for (int q = 0; q < n; q++) if (q == 0 || tmp[q] != tmp[q - 1]) UM[nrows++] = tmp[q];
        men[c] = (int)nrows;
    }
    free(tmp);
    R1 = malloc(sizeof(int) * (size_t)(npairs + 1)); R2 = malloc(sizeof(int) * (size_t)(npairs + 1));
    for (int c = 0; c < nU; c++)
        for (long long q = pst[c]; q < pen[c]; q++) {
            int lo_ = mst[c], hi_ = men[c] - 1;
            for (int s = 0; s < 2; s++) {
                int16_t m = s ? PM2[q] : PM1[q];
                int a = lo_, b = hi_;
                while (a < b) { int mid = (a + b) / 2; if (UM[mid] < m) a = mid + 1; else b = mid; }
                if (s) R2[q] = a; else R1[q] = a;
            }
        }
    /* FEAS, neighbours */
    for (int t = 0; t < nw; t++)
        for (int o = 0; o < 7; o++)
            for (int r = 0; r < 7; r++) {
                int ok = 0;
                for (int q = o; q <= o + r && q <= 6; q++) { int k = q < 7 - q ? q : 7 - q; if ((kall[t] >> k) & 1) ok = 1; }
                FEAS[t][o][r] = ok;
            }
    nwd = (nU + 63) / 64;
    nbr = calloc((size_t)nU * nwd, sizeof(u64));
    for (int i = 0; i < nU; i++)
        for (int j = i + 1; j < nU; j++) {
            if (!(__builtin_popcount(Fm[i] & cls[j]) & 1)) continue;
            if (idxof[cls[i] ^ cls[j]] < 0) continue;
            nbr[(size_t)i * nwd + (j >> 6)] |= 1ULL << (j & 63);
        }
    for (int d = 0; d < 7; d++) candl[d] = calloc(nwd, sizeof(u64));
    for (int i = 0; i < nU; i++) candl[0][i >> 6] |= 1ULL << (i & 63);
    memset(ocnt, 0, sizeof ocnt); memset(ucnt, 0, sizeof ucnt);
    memset(spanmark, 0, sizeof spanmark); spanmark[0] = 1; spanel[0] = 0; nspan = 1;
    fprintf(stderr, "vectors %d classes %d pairs %lld mask rows %lld witnesses %d (pairs %d)\n", nvec, nU, npairs, nrows, nw, NPW);
    /* phase 1: sampling over all witness pairs */
    killcnt = calloc(NPW + 1, sizeof(long long));
    alo = malloc(sizeof(int) * (NPW + 1)); ahi = malloc(sizeof(int) * (NPW + 1));
    acmn = malloc(sizeof(int) * (NPW + 1)); acmx = malloc(sizeof(int) * (NPW + 1));
    /* sampled kill sets for the greedy cover */
    int maxs = walks > 0 ? walks : 1;
    unsigned char *ks = calloc((size_t)maxs * NPW, 1);
    int nks = 0;
    for (int k = 0; k < walks; k++) {
        long long before = samp_killed;
        memset(killcnt, 0, 0);
        long long *snap = NULL; (void)snap;
        /* record kill set of this sample: diff of killcnt */
        static long long *prev = NULL;
        if (!prev) prev = calloc(NPW + 1, sizeof(long long));
        memcpy(prev, killcnt, sizeof(long long) * NPW);
        walk(0);
        if (samp_killed > before) {
            for (int p = 0; p < NPW; p++) if (killcnt[p] > prev[p]) ks[(size_t)nks * NPW + p] = 1;
            nks++;
        }
    }
    fprintf(stderr, "sampling: walks %d sampled leaves %lld killed %lld m1-dead %lld | est leaves %.4g est M1-valid %.4g\n",
            walks, samp_leaves, samp_killed, samp_m1, est_leaves / (walks ? walks : 1), est_m1 / (walks ? walks : 1));
    /* phase 2: greedy cover twice + frequency */
    Plist = malloc(sizeof(int) * (NPW + 1)); nP = 0;
    unsigned char *inP = calloc(NPW + 1, 1);
    for (int rnd = 0; rnd < 2; rnd++) {
        while (nP < pmax) {
            long long *cnt = calloc(NPW + 1, sizeof(long long));
            int any = 0;
            for (int s = 0; s < nks; s++) {
                int cov = 0;
                for (int p = 0; p < NPW; p++) if (inP[p] && ks[(size_t)s * NPW + p]) cov++;
                if (cov > rnd) continue;
                for (int p = 0; p < NPW; p++) if (!inP[p] && ks[(size_t)s * NPW + p]) { cnt[p]++; any = 1; }
            }
            if (!any) { free(cnt); break; }
            int best = -1;
            for (int p = 0; p < NPW; p++) if (!inP[p] && (best < 0 || cnt[p] > cnt[best])) best = p;
            free(cnt);
            if (best < 0) break;
            inP[best] = 1; Plist[nP++] = best;
        }
    }
    while (nP < pmax) {
        int best = -1;
        for (int p = 0; p < NPW; p++) if (!inP[p] && killcnt[p] > 0 && (best < 0 || killcnt[p] > killcnt[best])) best = p;
        if (best < 0) break;
        inP[best] = 1; Plist[nP++] = best;
    }
    target = malloc(sizeof(int) * (nP + 1));
    printf("PAIRS %d", nP);
    for (int p = 0; p < nP; p++) { target[p] = targetAll[Plist[p]]; printf(" %d,%d", TPa[Plist[p]], TPb[Plist[p]]); }
    printf("\n");
    /* phase 2b: tables */
    TMN = malloc((size_t)nrows * nP + 1); TMX = malloc((size_t)nrows * nP + 1);
    for (long long r = 0; r < nrows * nP; r++) { TMN[r] = 127; TMX[r] = -128; }
    for (int c = 0; c < nU; c++)
        for (long long q = pst[c]; q < pen[c]; q++) {
            const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
            for (int p = 0; p < nP; p++) {
                int P_ = Plist[p];
                int v = va[TPa[P_]] * va[TPb[P_]] + vb[TPa[P_]] * vb[TPb[P_]];
                if (v > 127 || v < -128) { puts("ERR int8"); return 1; }
                signed char *a1 = TMN + (size_t)R1[q] * nP + p, *b1 = TMX + (size_t)R1[q] * nP + p;
                signed char *a2 = TMN + (size_t)R2[q] * nP + p, *b2 = TMX + (size_t)R2[q] * nP + p;
                if (v < *a1) *a1 = (signed char)v;
                if (v > *b1) *b1 = (signed char)v;
                if (v < *a2) *a2 = (signed char)v;
                if (v > *b2) *b2 = (signed char)v;
            }
        }
    lo = malloc(sizeof(int) * (nP + 1)); hi = malloc(sizeof(int) * (nP + 1));
    cmn = malloc(sizeof(int) * (nP + 1)); cmx = malloc(sizeof(int) * (nP + 1));
    if (!strcmp(argv[2], "est")) {
        est_nodes = est_leaves = est_surv = 0;
        int w2 = walks;
        for (int k = 0; k < w2; k++) walk(1);
        printf("KNUTH walks %d est_nodes %.4g est_leaves %.4g est_survivors %.4g\n", w2, est_nodes / w2, est_leaves / w2,
               est_surv / w2);
        return 0;
    }
    fout = strcmp(argv[5], "-") ? fopen(argv[5], "w") : NULL;
    for (int s = 1; s < nspan; s++) spanmark[spanel[s]] = 0;
    memset(spanmark, 0, sizeof spanmark); spanmark[0] = 1; spanel[0] = 0; nspan = 1;
    memset(ocnt, 0, sizeof ocnt); memset(ucnt, 0, sizeof ucnt);
    rec(0);
    printf("DONE nodes %lld leaves %lld m1_dead %lld budget_dead %lld survivors %lld\n", nodes, leaves, m1dead, bdead, surv);
    if (fout) fclose(fout);
    return 0;
}
