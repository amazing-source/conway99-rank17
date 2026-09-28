/* skel_w.c -- skeleton enumeration with witness budgets at the leaves (structures, 28/09 night).
   Skeleton = 6-set {g_2..g_7} of usable classes of K/2K (search order = index order of the exported classes), pairwise odd
   forms, all g_c ^ g_d usable, independent (C1's rules), first element an H-orbit minimum (C1's H-reduction), and
   (new, S7 base-line reduction) g_2 = the first element is the smallest of all 21 cell classes in index order: for any
   frame, re-basing at an end of the cell holding the smallest class puts it in S_1, and the H-move to its orbit minimum
   keeps it smallest (index order sorts by (orbit min, class)).
   Witness o-rules: o_t = #{d : (kappa_t, g_d) odd}, k_t = min(o_t, 7 - o_t) must lie in the allowed set of witness t.
   Unit rule: every unit functional is odd on some g_d.
   Leaf check (M1-aware budgets): cell (1,d) -> class g_d with M1 pattern all ones; cell (c,d), 2 <= c < d -> class
   g_c ^ g_d with pattern e_c (over the six g's).  A class carries, per a-mask m (functional of a = (w+ +- w-)/2 mod 2 over
   the valid pairs), int8 minima / maxima of the contributions a_t(w+)a_t'(w+) + a_t(w-)a_t'(w-) for nP witness pairs.
   For each cell, the masks whose pattern ((m, g_e) mod 2)_e equals the cell pattern are allowed (none -> dead: M1);
   budget: sum of cell minima <= target <= sum of cell maxima for every exported pair.
   Input (binary, little endian): int32 nU, nw, nunits, nP; int32 target[nP]; int8 kallowed[nw] (bit k set = k allowed);
     per class: int32 cls, int32 Fm, int8 orbmin, uint64 wpar, uint32 upar, int32 nm, int16 masks[nm], int8 mn[nm*nP],
     int8 mx[nm*nP].
   Usage: skel_w.exe input survivors_out nshards shard [knuth_walks seed]
   Output (stdout): nodes leaves m1_dead budget_dead survivors (+ Knuth estimates in walk mode).                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
static int nU, nw, nun, nP, nwd;
static int *target, *cls, *Fm, *nm;
static signed char *orbmin, *kall;
static u64 *wpar;
static uint32_t *upar;
static int16_t **masks;
static signed char **mn, **mx;
static int idxof[1 << 15];
static u64 *nbr;                      /* nU * nwd */
static u64 **candl;                   /* per depth */
static int FEAS[64][7][7];            /* per witness t (<64): (o so far, r left) feasible */
static long long nodes = 0, leaves = 0, m1dead = 0, bdead = 0, surv = 0;
static FILE *fout;
static int chosen[6];
static int nshards = 1, shard = 0;
static int *lo, *hi, *cmn, *cmx;

static int pc(u64 x) { return __builtin_popcountll(x); }

static int leaf_check(void) {
    int G[6];
    for (int e = 0; e < 6; e++) G[e] = cls[chosen[e]];
    for (int p = 0; p < nP; p++) { lo[p] = 0; hi[p] = 0; }
    for (int c = -1; c < 6; c++)
        for (int d = c + 1; d < 6; d++) {
            int C = (c < 0) ? G[d] : (G[c] ^ G[d]);
            int i = idxof[C];
            if (i < 0) return 0;
            int want = (c < 0) ? 63 : (1 << c);
            int any = 0;
            for (int p = 0; p < nP; p++) { cmn[p] = 1 << 20; cmx[p] = -(1 << 20); }
            for (int q = 0; q < nm[i]; q++) {
                int m = masks[i][q], pat = 0;
                for (int e = 0; e < 6; e++) pat |= (__builtin_popcount(m & G[e]) & 1) << e;
                if (pat != want) continue;
                any = 1;
                const signed char *a = mn[i] + (size_t)q * nP, *b = mx[i] + (size_t)q * nP;
                for (int p = 0; p < nP; p++) { if (a[p] < cmn[p]) cmn[p] = a[p]; if (b[p] > cmx[p]) cmx[p] = b[p]; }
            }
            if (!any) { m1dead++; return 0; }
            for (int p = 0; p < nP; p++) { lo[p] += cmn[p]; hi[p] += cmx[p]; }
        }
    for (int p = 0; p < nP; p++)
        if (lo[p] > target[p] || hi[p] < target[p]) { bdead++; return 0; }
    return 1;
}

/* span marks for independence */
static unsigned char spanmark[1 << 15];
static int spanel[64], nspan = 1;

static int ocnt[7][64], ucnt[7][32];

/* returns number of valid children at this node when counting (Knuth), else recurses */
static void rec(int depth) {
    nodes++;
    u64 *cand = candl[depth];
    int need = 5 - depth;
    for (int w = 0; w < nwd; w++) {
        u64 x = cand[w];
        while (x) {
            int bit = __builtin_ctzll(x); x &= x - 1;
            int i = 64 * w + bit;
            int g = cls[i];
            if (spanmark[g]) continue;
            if (depth == 0) {
                if (!orbmin[i]) continue;
                if (nshards > 1 && (i % nshards) != shard) continue;
            } else {
                int bad = 0;
                for (int j = 1; j < depth && !bad; j++) if (idxof[cls[chosen[j]] ^ g] < chosen[0]) bad = 1;
                if (idxof[cls[chosen[0]] ^ g] < chosen[0]) bad = 1;
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
                    if (fout) { for (int e = 0; e < 6; e++) fprintf(fout, "%d ", cls[chosen[e]]); fprintf(fout, "\n"); fflush(fout); }
                }
                continue;
            }
            u64 *c2 = candl[depth + 1], *nb = nbr + (size_t)i * nwd;
            int cnt = 0;
            for (int k = 0; k < nwd; k++) { c2[k] = cand[k] & nb[k]; cnt += pc(c2[k]); }
            if (cnt < need) continue;
            /* extend span */
            int old = nspan;
            for (int s = 0; s < old; s++) { int y = spanel[s] ^ g; spanel[nspan++] = y; spanmark[y] = 1; }
            rec(depth + 1);
            for (int s = old; s < nspan; s++) spanmark[spanel[s]] = 0;
            nspan = old;
        }
    }
}

/* Knuth random probe: at each depth collect valid children, pick one uniformly, multiply */
static unsigned long long rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static double est_nodes, est_leaves, est_surv;
static int kids[1 << 15];
static FILE *sampf = NULL;

static void walk(void) {
    double mult = 1.0;
    nspan = 1; spanel[0] = 0; memset(spanmark, 0, sizeof spanmark); spanmark[0] = 1;
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
                for (int t = 0; t < nw; t++) {
                    int o = ocnt[depth][t] + (int)((wpar[i] >> t) & 1);
                    if (!FEAS[t][o][r]) { bad = 1; break; }
                }
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
            /* evaluate every leaf child exactly (cheap relative to the walk) */
            int ns = 0, first = -1;
            for (int q = 0; q < nk; q++) {
                int qq = (int)((q + (sampf ? xr() % nk : 0)) % nk);
                chosen[5] = kids[qq];
                if (leaf_check()) { ns++; if (first < 0) first = kids[qq]; }
            }
            est_leaves += mult * nk; est_surv += mult * ns;
            if (sampf && first >= 0) {
                chosen[5] = first;
                for (int e = 0; e < 6; e++) fprintf(sampf, "%d ", cls[chosen[e]]);
                fputc(10, sampf);
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

int main(int argc, char **argv) {
    if (argc < 5) { puts("usage: skel_w.exe input survivors nshards shard [walks seed]"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { puts("ERR open"); return 1; }
    int hdr[4];
    if (fread(hdr, 4, 4, f) != 4) return 1;
    nU = hdr[0]; nw = hdr[1]; nun = hdr[2]; nP = hdr[3];
    if (nw > 64 || nun > 32) { puts("ERR sizes"); return 1; }
    target = malloc(sizeof(int) * (nP + 1));
    if (nP && fread(target, 4, nP, f) != (size_t)nP) return 1;
    kall = malloc(nw + 1);
    if (nw && fread(kall, 1, nw, f) != (size_t)nw) return 1;
    cls = malloc(sizeof(int) * nU); Fm = malloc(sizeof(int) * nU); orbmin = malloc(nU); wpar = malloc(8 * (size_t)nU);
    upar = malloc(4 * (size_t)nU); nm = malloc(sizeof(int) * nU);
    masks = malloc(sizeof(int16_t *) * nU); mn = malloc(sizeof(signed char *) * nU); mx = malloc(sizeof(signed char *) * nU);
    for (int i = 0; i < (1 << 15); i++) idxof[i] = -1;
    for (int i = 0; i < nU; i++) {
        if (fread(&cls[i], 4, 1, f) != 1 || fread(&Fm[i], 4, 1, f) != 1 || fread(&orbmin[i], 1, 1, f) != 1 ||
            fread(&wpar[i], 8, 1, f) != 1 || fread(&upar[i], 4, 1, f) != 1 || fread(&nm[i], 4, 1, f) != 1) { puts("ERR class"); return 1; }
        masks[i] = malloc(2 * (size_t)(nm[i] + 1));
        mn[i] = malloc((size_t)nm[i] * nP + 1); mx[i] = malloc((size_t)nm[i] * nP + 1);
        if (nm[i] && (fread(masks[i], 2, nm[i], f) != (size_t)nm[i] || fread(mn[i], 1, (size_t)nm[i] * nP, f) != (size_t)nm[i] * nP ||
                      fread(mx[i], 1, (size_t)nm[i] * nP, f) != (size_t)nm[i] * nP)) { puts("ERR tables"); return 1; }
        idxof[cls[i]] = i;
    }
    fclose(f);
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
    candl = malloc(sizeof(u64 *) * 7);
    for (int d = 0; d < 7; d++) candl[d] = calloc(nwd, sizeof(u64));
    for (int i = 0; i < nU; i++) candl[0][i >> 6] |= 1ULL << (i & 63);
    lo = malloc(sizeof(int) * (nP + 1)); hi = malloc(sizeof(int) * (nP + 1));
    cmn = malloc(sizeof(int) * (nP + 1)); cmx = malloc(sizeof(int) * (nP + 1));
    memset(ocnt, 0, sizeof ocnt); memset(ucnt, 0, sizeof ucnt);
    nshards = atoi(argv[3]); shard = atoi(argv[4]);
    if (argc > 5) {
        int walks = atoi(argv[5]);
        if (argc > 6) rng ^= (unsigned long long)atoll(argv[6]) * 0x9E3779B97F4A7C15ULL;
        if (argc > 7) sampf = fopen(argv[7], "w");
        for (int k = 0; k < walks; k++) walk();
        printf("KNUTH walks %d est_nodes %.4g est_leaves %.4g est_survivors %.4g budget_dead %lld m1_dead %lld\n", walks,
               est_nodes / walks, est_leaves / walks, est_surv / walks, bdead, m1dead);
        if (sampf) fclose(sampf);
        return 0;
    }
    fout = strcmp(argv[2], "-") ? fopen(argv[2], "w") : NULL;
    spanel[0] = 0; spanmark[0] = 1; nspan = 1;
    rec(0);
    printf("DONE nodes %lld leaves %lld m1_dead %lld budget_dead %lld survivors %lld\n", nodes, leaves, m1dead, bdead, surv);
    if (fout) fclose(fout);
    return 0;
}
