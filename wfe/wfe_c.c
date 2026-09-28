/* wfe_c.c -- per-class witness-filtered frame engine, all levels in C (structures, 28/09 night).
   Skeleton level (skel_w3.c): C1's skeleton rules (usable classes, pairwise odd forms, usable sums, independence, units),
   witness o-rules (k_t = min(o_t, 7 - o_t) allowed), S7 base-line reduction (g_2 = smallest of the 21 cell classes),
   lexmin stabilizer pruning (t_{j+1} minimal in its orbit under the transvections fixing t_1..t_j), table pre-check
   (M1-aware min/max over a sampled pair subset P).
   Leaf level (exact, per skeleton), for every cell (c,d) (line 0 = base line 1, line e+1 = g_{e+2}):
     pairs {w, w'} of the cell class with (w, w') in {0, +-4} (filtered shell: per-vector witness filter, basic L1);
     M1: a = (w + w')/2 or (w - w')/2 has pattern ((a, g_e) mod 2)_e = eps(c) (eps(line 0) = all ones, eps(e+1) = e_e);
     refined L1 (C1): w bad if some v = (r - w)/2 (r a root of V, r = w mod 2V) has pattern in {0, eps_c, eps_d,
       eps_c + eps_d} up to complement;
     exact-k values: n=6 witness with k=3 -> no |a| = 2; n=6 with k=2 -> not both labels |a| = 2; n=7 with k=3 -> no
       |a| = 3 (the other value rules are k-independent and applied to the shell);
     budgets: for every witness pair (t,t'), sum over cells of a_t a_t' (both labels) = 28 (kappa_t, kappa_t') has to lie
       between the sums of cell minima and maxima.
   Vector level: signed ordered entries (w+, w-) with a = (w+ + w-)/2 of the M1 pattern, sign normalisation (a_12, b_12,
     b_1d first non-zero K-coordinate positive), then the sbdfs_w search in K coordinates: pairwise rule
     ((w_X, w_Y) + r_X.r_Y)/4 in {0, +-1} (lazy forward checking), star closure + look-ahead (1 or 2 open cells), PSD
     residual 28 I - sum proj proj^T (proj = L^T w, GK = L L^T), witness budgets (running sums, static suffix bounds).
   Input (binary LE): int32 nvec, nw, nU, nun; int32 GK[225]; int32 target[nw(nw+1)/2]; int8 kall[nw]; int8 wn[nw];
     per vector: int8 x[15], int8 val[nw], int32 cls_index; per class: int32 cls, int32 Fm, int8 orbmin, uint64 wpar,
     uint32 upar; int32 nroots, nroots x (int32 rb, int32 rf); per vector: int32 cnt, int16 l1mask[cnt].
   Usage: wfe_c.exe input est|run nshards shard outprefix walks seed [pmax] [nodecap]
   Env: STABMAX (default 5).                                                                                             */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef uint64_t u64;
#define NC 21
static int nvec, nw, nU, nun, NPW, nP, nwd, pmax = 40;
static long long nodecap = 20000000;
static int GK[15][15];
static double Lt[15][15];
static int *targetAll, *target, *Plist;
static signed char *kall, *wn, *X, *VAL;
static int *vcls, *vstart, *vend;
static u64 *M26, *M37;
static int *l1off; static int16_t *l1m;
static int *cls, *Fm; static signed char *orbmin; static u64 *wpar; static uint32_t *upar;
static int idxof[1 << 15];
static int chosen[6];
static int nroots = 0, *RB, *RF, STABMAX = 5;
static int *uf[7];
static long long stab_pruned = 0;
static int *pst, *pen, *PI, *PJ; static int16_t *PM1, *PM2; static long long npairs = 0;
static int *mst, *men; static int16_t *UM; static int *R1, *R2;
static long long *RPoff; static int *RPq; static signed char *RPo;
static signed char *TMN, *TMX;
static int *TPa, *TPb;
static u64 *nbr; static u64 *candl[7];
static int FEAS[64][7][7];
static long long nodes = 0, leaves = 0, m1dead = 0, bdead = 0, pass_tab = 0;
static long long ex_ac = 0; static int pair_ac = 0;
static long long ex_empty = 0, ex_budget = 0, sv_excl = 0, sv_limit = 0, sv_found = 0, sv_nodes_tot = 0, sv_root = 0;
static int nshards = 1, shard = 0;
static FILE *flim, *ffound;
static int ocnt[7][64], ucnt[7][32];
static unsigned char spanmark[1 << 15];
static int spanel[64], nspan = 1;
static long long *killcnt; static int *lo, *hi, *cmn, *cmx;
static unsigned long long rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int pc(u64 x) { return __builtin_popcountll(x); }
static int patof(int m, const int *G) { int p = 0; for (int e = 0; e < 6; e++) p |= (__builtin_popcount(m & G[e]) & 1) << e; return p; }
static inline int dotxg(int u, int v);
static int ufind(int *p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

static void build_orbits(int depth) {
    int *p = uf[depth];
    for (int i = 0; i < nU; i++) p[i] = i;
    for (int q = 0; q < nroots; q++) {
        int ok = 1;
        for (int j = 0; j < depth && ok; j++) if (__builtin_popcount(cls[chosen[j]] & RF[q]) & 1) ok = 0;
        if (!ok) continue;
        for (int i = 0; i < nU; i++) {
            int x = cls[i];
            if (!(__builtin_popcount(x & RF[q]) & 1)) continue;
            int j = idxof[x ^ RB[q]];
            if (j < 0) continue;
            int a = ufind(p, i), b = ufind(p, j);
            if (a < b) p[b] = a; else if (b < a) p[a] = b;
        }
    }
}

/* ------------------------------------------------------------------ table pre-check (P subset) */
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

/* ------------------------------------------------------------------ exact leaf check */
typedef struct { int q; signed char o1, o2; } CP;
static CP *cellp[NC]; static int ncellp[NC], cellcap[NC];
static int celli[NC], cellc[NC], celld[NC];          /* class index, lines 0..6 */
static int LG[6];                                      /* the six g's (class bits) */

static int eps_of(int line) { return line == 0 ? 63 : 1 << (line - 1); }

static int exact_check(void) {
    for (int e = 0; e < 6; e++) LG[e] = cls[chosen[e]];
    /* witness k's */
    u64 K36 = 0, K26 = 0, K37 = 0;
    for (int t = 0; t < nw; t++) {
        int o = 0;
        for (int e = 0; e < 6; e++) o += (int)((wpar[chosen[e]] >> t) & 1);
        int k = o < 7 - o ? o : 7 - o;
        if (wn[t] == 6 && k == 3) K36 |= 1ULL << t;
        if (wn[t] == 6 && k == 2) K26 |= 1ULL << t;
        if (wn[t] == 7 && k == 3) K37 |= 1ULL << t;
    }
    int x = 0;
    for (int c = 0; c < 7; c++)
        for (int d = c + 1; d < 7; d++) {
            int C = (c == 0) ? LG[d - 1] : (LG[c - 1] ^ LG[d - 1]);
            int i = idxof[C];
            celli[x] = i; cellc[x] = c; celld[x] = d;
            int want = eps_of(c), ec = eps_of(c), ed = eps_of(d);
            int bad[8] = {0, ec, ed, ec ^ ed, 63, 63 ^ ec, 63 ^ ed, 63 ^ ec ^ ed};
            int n = 0;
            for (int rw = mst[i]; rw < men[i]; rw++) {
              if (patof(UM[rw], LG) != want) continue;
              for (long long z_ = RPoff[rw]; z_ < RPoff[rw + 1]; z_++) {
                int q = RPq[z_];
                int o1 = RPo[z_] == 1, o2 = RPo[z_] == 2;
                int a = PI[q], b = PJ[q];
                if ((M26[a] | M26[b]) & K36) continue;
                if (M26[a] & M26[b] & K26) continue;
                if ((M37[a] | M37[b]) & K37) continue;
                int l1bad = 0;
                for (int s = 0; s < 2 && !l1bad; s++) {
                    int v = s ? b : a;
                    for (int j = l1off[v]; j < l1off[v + 1] && !l1bad; j++) {
                        int p = patof(l1m[j], LG);
                        for (int z = 0; z < 8; z++) if (p == bad[z]) { l1bad = 1; break; }
                    }
                }
                if (l1bad) continue;
                if (n == cellcap[x]) { cellcap[x] = cellcap[x] ? 2 * cellcap[x] : 256; cellp[x] = realloc(cellp[x], sizeof(CP) * cellcap[x]); }
                cellp[x][n].q = (int)q; cellp[x][n].o1 = (signed char)o1; cellp[x][n].o2 = (signed char)o2; n++;
              }
            }
            ncellp[x] = n;
            if (!n) { ex_empty++; return 0; }
            x++;
        }
    /* budgets over all witness pairs, fixpoint propagation on the pairs */
    static signed char *PCB; static size_t pcbcap = 0; static char *alive; static size_t alcap = 0;
    int off[NC + 1]; off[0] = 0;
    for (x = 0; x < NC; x++) off[x + 1] = off[x] + ncellp[x];
    int tot = off[NC];
    if ((size_t)tot * NPW > pcbcap) { pcbcap = (size_t)tot * NPW * 2; PCB = realloc(PCB, pcbcap); }
    if ((size_t)tot > alcap) { alcap = (size_t)tot * 2; alive = realloc(alive, alcap); }
    for (x = 0; x < NC; x++)
        for (int k = 0; k < ncellp[x]; k++) {
            int q = cellp[x][k].q;
            const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
            signed char *out = PCB + (size_t)(off[x] + k) * NPW;
            for (int p = 0; p < NPW; p++) out[p] = (signed char)(va[TPa[p]] * va[TPb[p]] + vb[TPa[p]] * vb[TPb[p]]);
            alive[off[x] + k] = 1;
        }
    static int *pmn, *pmx; static int pmcap = 0;
    if (pmcap < NC * NPW) { pmcap = NC * NPW; pmn = realloc(pmn, sizeof(int) * pmcap); pmx = realloc(pmx, sizeof(int) * pmcap); }
    for (int pass = 0; pass < 40; pass++) {
        for (int p = 0; p < NPW; p++) { alo[p] = 0; ahi[p] = 0; }
        for (x = 0; x < NC; x++) {
            int *mn = pmn + x * NPW, *mx = pmx + x * NPW, any = 0;
            for (int p = 0; p < NPW; p++) { mn[p] = 1 << 20; mx[p] = -(1 << 20); }
            for (int k = off[x]; k < off[x + 1]; k++) {
                if (!alive[k]) continue;
                any = 1;
                const signed char *c = PCB + (size_t)k * NPW;
                for (int p = 0; p < NPW; p++) { if (c[p] < mn[p]) mn[p] = c[p]; if (c[p] > mx[p]) mx[p] = c[p]; }
            }
            if (!any) { ex_budget++; return 0; }
            for (int p = 0; p < NPW; p++) { alo[p] += mn[p]; ahi[p] += mx[p]; }
        }
        for (int p = 0; p < NPW; p++) if (alo[p] > targetAll[p] || ahi[p] < targetAll[p]) { ex_budget++; return 0; }
        int changed = 0;
        for (x = 0; x < NC; x++) {
            int *mn = pmn + x * NPW, *mx = pmx + x * NPW;
            for (int k = off[x]; k < off[x + 1]; k++) {
                if (!alive[k]) continue;
                const signed char *c = PCB + (size_t)k * NPW;
                for (int p = 0; p < NPW; p++)
                    if (c[p] + alo[p] - mn[p] > targetAll[p] || c[p] + ahi[p] - mx[p] < targetAll[p]) { alive[k] = 0; changed = 1; break; }
            }
        }
        if (!changed) break;
    }
    for (x = 0; x < NC; x++) {
        int n = 0;
        for (int k = 0; k < ncellp[x]; k++) if (alive[off[x] + k]) cellp[x][n++] = cellp[x][k];
        ncellp[x] = n;
        if (!n) { ex_budget++; return 0; }
    }
    if (pair_ac) {
        /* pair-level arc consistency: disjoint cells need all four |(w, w')| in {0, 4}; cells sharing a line need
           |(w, w')| <= 5 (w, w' vectors of the two pairs) */
        int changed = 1, rounds = 0;
        while (changed && rounds++ < 20) {
            changed = 0;
            for (int X_ = 0; X_ < NC; X_++) {
                int n = 0;
                for (int k = 0; k < ncellp[X_]; k++) {
                    int q = cellp[X_][k].q, a1 = PI[q], a2 = PJ[q], ok = 1;
                    for (int Y = 0; Y < NC && ok; Y++) {
                        if (Y == X_) continue;
                        int shared = (cellc[X_] == cellc[Y] || cellc[X_] == celld[Y] || celld[X_] == cellc[Y] || celld[X_] == celld[Y]);
                        int sup = 0;
                        for (int m = 0; m < ncellp[Y] && !sup; m++) {
                            int r = cellp[Y][m].q, b1 = PI[r], b2 = PJ[r];
                            int i1 = dotxg(a1, b1), i2 = dotxg(a1, b2), i3 = dotxg(a2, b1), i4 = dotxg(a2, b2);
                            if (shared) {
                                if (abs(i1) <= 5 && abs(i2) <= 5 && abs(i3) <= 5 && abs(i4) <= 5) sup = 1;
                            } else {
                                if ((i1 == 0 || abs(i1) == 4) && (i2 == 0 || abs(i2) == 4) && (i3 == 0 || abs(i3) == 4) && (i4 == 0 || abs(i4) == 4)) sup = 1;
                            }
                        }
                        if (!sup) ok = 0;
                    }
                    if (ok) cellp[X_][n++] = cellp[X_][k]; else changed = 1;
                }
                ncellp[X_] = n;
                if (!n) { ex_ac++; return 0; }
            }
        }
    }
    return 1;
}

/* ------------------------------------------------------------------ vector-level solver (K coordinates) */
/* entry = (w+, w-) = (sa X[va], sb X[vb]); per shell vector: GWV = GK x (int), LWV = L^T x (double) */
typedef struct { int va, vb, q; signed char sa, sb; u64 ka[2], kb[2]; } ENT;
static int *GWV; static double *LWV;
static ENT *ent[NC]; static int nent[NC], entcap[NC];
static signed char *econ[NC];
static int ncons, *ctarget, *cpair, *csmin, *csmax;
static int ord[NC], CCo[NC], CDo[NC], lastk[7], rrt[NC][NC][2][2];
static u64 *lev[NC + 1][NC]; static int nwo[NC];
static int schosen[NC];
static long long snodes, dprof[NC + 1]; static int crun[NC + 1][4096];
static int *htab; static int hcap = 0;
static int dynb = 0, dlo[4096], dhi[4096], dmn[4096], dmx[4096];

static void packv(const int *v, u64 *out) {
    unsigned char b[16]; memset(b, 0, 16);
    for (int i = 0; i < 15; i++) b[i] = (unsigned char)(v[i] + 128);
    memcpy(out, b, 16);
}
static u64 hash2(const u64 *k) { u64 h = k[0] * 0x9E3779B97F4A7C15ULL ^ (k[1] + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL; return h ^ (h >> 29); }
static int fpos(const int *v) { for (int i = 0; i < 15; i++) { if (v[i] > 0) return 1; if (v[i] < 0) return 0; } return 1; }

static void add_entry(int x, int va, int sa, int vb, int sb, int q) {
    int c = cellc[x], d = celld[x];
    const signed char *A = X + (size_t)va * 15, *B = X + (size_t)vb * 15;
    int s2a[15], s2b[15], a[15], b[15];
    for (int i = 0; i < 15; i++) { s2a[i] = sa * A[i] + sb * B[i]; s2b[i] = sa * A[i] - sb * B[i]; a[i] = s2a[i] / 2; b[i] = s2b[i] / 2; }
    if (c == 0 && d == 1) { if (!fpos(a) || !fpos(b)) return; }
    else if (c == 0) { if (!fpos(b)) return; }
    if (nent[x] == entcap[x]) { entcap[x] = entcap[x] ? 2 * entcap[x] : 512; ent[x] = realloc(ent[x], sizeof(ENT) * entcap[x]); }
    ENT *E = &ent[x][nent[x]++];
    E->va = va; E->vb = vb; E->sa = (signed char)sa; E->sb = (signed char)sb; E->q = q;
    packv(s2a, E->ka); packv(s2b, E->kb);
}

static void starsum(int v, int k, int *s) {
    memset(s, 0, 15 * sizeof(int));
    for (int x = 0; x <= k; x++) {
        const ENT *E = &ent[ord[x]][schosen[x]];
        const signed char *A = X + (size_t)E->va * 15, *B = X + (size_t)E->vb * 15;
        if (CCo[x] == v) for (int i = 0; i < 15; i++) s[i] += E->sa * A[i] + E->sb * B[i];
        else if (CDo[x] == v) for (int i = 0; i < 15; i++) s[i] += E->sa * A[i] - E->sb * B[i];
    }
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

static int anybit(const u64 *b, int n) { for (int i = 0; i < n; i++) if (b[i]) return 1; return 0; }

static inline int dotxg(int u, int v) {          /* X[u] . GK X[v] */
    const signed char *a = X + (size_t)u * 15; const int *g = GWV + (size_t)v * 15; int s = 0;
    for (int i = 0; i < 15; i++) s += a[i] * g[i];
    return s;
}

static int compat(int Xp, int t, int Yp, int q) {
    const ENT *A = &ent[ord[Xp]][t], *B = &ent[ord[Yp]][q];
    int va[2] = {A->va, A->vb}, sa[2] = {A->sa, A->sb}, vb[2] = {B->va, B->vb}, sb[2] = {B->sa, B->sb};
    for (int s = 0; s < 2; s++)
        for (int u = 0; u < 2; u++) {
            int v = rrt[Xp][Yp][s][u] + sa[s] * sb[u] * dotxg(va[s], vb[u]);
            if (!(v == 0 || v == 4 || v == -4)) return 0;
        }
    return 1;
}

static int solve_rec(int k, double R[15][15]);

static int try_entry(int k, int t, double R[15][15]) {
    int x = ord[k];
    for (int i = 0; i < ncons; i++) {
        int v = crun[k][i] + econ[x][(size_t)t * ncons + i];
        if (v + csmin[(k + 1) * ncons + i] > ctarget[i] || v + csmax[(k + 1) * ncons + i] < ctarget[i]) return 0;
        crun[k + 1][i] = v;
    }
    const ENT *E = &ent[x][t];
    const double *pa = LWV + (size_t)E->va * 15, *pb = LWV + (size_t)E->vb * 15;
    double R2[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pa[i] * pa[j] - pb[i] * pb[j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    schosen[k] = t;
    int s[15];
    for (int v = 0; v < 7; v++)
        if (lastk[v] == k) { starsum(v, k, s); for (int i = 0; i < 15; i++) if (s[i]) return 0; }
    for (int Y = k + 1; Y < NC; Y++) {
        u64 *src = lev[k][Y], *dst = lev[k + 1][Y], acc = 0;
        for (int w = 0; w < nwo[Y]; w++) {
            u64 xx = src[w], keep = 0;
            while (xx) {
                int bit = __builtin_ctzll(xx); xx &= xx - 1;
                if (compat(k, t, Y, 64 * w + bit)) keep |= 1ULL << bit;
            }
            dst[w] = keep; acc |= keep;
        }
        if (!acc) return 0;
    }
    if (dynb && k < dynb) {
        for (int i = 0; i < ncons; i++) { dlo[i] = crun[k + 1][i]; dhi[i] = crun[k + 1][i]; }
        for (int Y = k + 1; Y < NC; Y++) {
            int y = ord[Y];
            for (int i = 0; i < ncons; i++) { dmn[i] = 1 << 20; dmx[i] = -(1 << 20); }
            for (int w = 0; w < nwo[Y]; w++) {
                u64 xx = lev[k + 1][Y][w];
                while (xx) {
                    int bit = __builtin_ctzll(xx); xx &= xx - 1;
                    const signed char *ec = econ[y] + (size_t)(64 * w + bit) * ncons;
                    for (int i = 0; i < ncons; i++) { if (ec[i] < dmn[i]) dmn[i] = ec[i]; if (ec[i] > dmx[i]) dmx[i] = ec[i]; }
                }
            }
            for (int i = 0; i < ncons; i++) { dlo[i] += dmn[i]; dhi[i] += dmx[i]; }
        }
        for (int i = 0; i < ncons; i++) if (dlo[i] > ctarget[i] || dhi[i] < ctarget[i]) return 0;
    }
    for (int v = 0; v < 7; v++) {
        int op[3], cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CCo[Y] == v || CDo[Y] == v) { if (cnt < 3) op[cnt] = Y; cnt++; }
        if (cnt == 1) {
            int Y = op[0];
            starsum(v, k, s);
            u64 *b = lev[k + 1][Y];
            for (int w = 0; w < nwo[Y]; w++) {
                u64 xx = b[w];
                while (xx) {
                    int bit = __builtin_ctzll(xx); xx &= xx - 1;
                    const ENT *U = &ent[ord[Y]][64 * w + bit];
                    const signed char *A = X + (size_t)U->va * 15, *B = X + (size_t)U->vb * 15;
                    for (int i = 0; i < 15; i++) {
                        int c = (CCo[Y] == v) ? U->sa * A[i] + U->sb * B[i] : U->sa * A[i] - U->sb * B[i];
                        if (c != -s[i]) { b[w] &= ~(1ULL << bit); break; }
                    }
                }
            }
            if (!anybit(b, nwo[Y])) return 0;
        } else if (cnt == 2) {
            starsum(v, k, s);
            for (int side = 0; side < 2; side++) {
                int Ac = op[side], Bc = op[1 - side];
                u64 *bB = lev[k + 1][Bc], *bA = lev[k + 1][Ac];
                int nB = 0;
                for (int w = 0; w < nwo[Bc]; w++) nB += pc(bB[w]);
                int hs = 1;
                while (hs < 2 * nB + 2) hs <<= 1;
                if (hs > hcap) { free(htab); htab = malloc(sizeof(int) * (size_t)hs); hcap = hs; }
                memset(htab, 0, sizeof(int) * (size_t)hs);
                for (int w = 0; w < nwo[Bc]; w++) {
                    u64 xx = bB[w];
                    while (xx) {
                        int bit = __builtin_ctzll(xx); xx &= xx - 1;
                        int q = 64 * w + bit;
                        const ENT *U = &ent[ord[Bc]][q];
                        const u64 *kk = (CCo[Bc] == v) ? U->ka : U->kb;
                        int pos = (int)(hash2(kk) & (u64)(hs - 1));
                        while (htab[pos]) pos = (pos + 1) & (hs - 1);
                        htab[pos] = q + 1;
                    }
                }
                for (int w = 0; w < nwo[Ac]; w++) {
                    u64 xx = bA[w];
                    while (xx) {
                        int bit = __builtin_ctzll(xx); xx &= xx - 1;
                        const ENT *U = &ent[ord[Ac]][64 * w + bit];
                        const signed char *A = X + (size_t)U->va * 15, *B = X + (size_t)U->vb * 15;
                        int need[15], ok = 1;
                        for (int i = 0; i < 15; i++) {
                            int c = (CCo[Ac] == v) ? U->sa * A[i] + U->sb * B[i] : U->sa * A[i] - U->sb * B[i];
                            need[i] = -s[i] - c;
                            if (need[i] < -127 || need[i] > 127) ok = 0;
                        }
                        int found = 0;
                        if (ok) {
                            u64 key[2]; packv(need, key);
                            int pos = (int)(hash2(key) & (u64)(hs - 1));
                            while (htab[pos]) {
                                const ENT *Q = &ent[ord[Bc]][htab[pos] - 1];
                                const u64 *kk = (CCo[Bc] == v) ? Q->ka : Q->kb;
                                if (kk[0] == key[0] && kk[1] == key[1]) { found = 1; break; }
                                pos = (pos + 1) & (hs - 1);
                            }
                        }
                        if (!found) bA[w] &= ~(1ULL << bit);
                    }
                }
                if (!anybit(bA, nwo[Ac])) return 0;
            }
        }
    }
    return solve_rec(k + 1, R2);
}

static int solve_rec(int k, double R[15][15]) {
    snodes++;
    dprof[k]++;
    if (snodes > nodecap) return 2;
    if (k == NC) {
        for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0;
        return 1;
    }
    u64 *b = lev[k][k];
    for (int w = 0; w < nwo[k]; w++) {
        u64 xx = b[w];
        while (xx) {
            int bit = __builtin_ctzll(xx); xx &= xx - 1;
            int r = try_entry(k, 64 * w + bit, R);
            if (r) return r;
        }
    }
    return 0;
}

static int cmpsz(const void *a, const void *b) { return nent[*(const int *)a] - nent[*(const int *)b]; }

static int vector_level(void) {
    for (int x = 0; x < NC; x++) {
        nent[x] = 0;
        for (int k = 0; k < ncellp[x]; k++) {
            int q = cellp[x][k].q, va = PI[q], vb = PJ[q];
            for (int o = 0; o < 2; o++) {
                if (!(o ? cellp[x][k].o2 : cellp[x][k].o1)) continue;
                int s = o ? -1 : 1;                  /* a = (w_a + s w_b)/2 */
                add_entry(x, va, 1, vb, s, q); add_entry(x, vb, s, va, 1, q);
                add_entry(x, va, -1, vb, -s, q); add_entry(x, vb, -s, va, -1, q);
            }
        }
        if (!nent[x]) { sv_root++; return 0; }
    }
    for (int k = 0; k < NC; k++) ord[k] = k;
    if (!getenv("ORDERLEX")) qsort(ord, NC, sizeof(int), cmpsz);
    for (int k = 0; k < NC; k++) { CCo[k] = cellc[ord[k]]; CDo[k] = celld[ord[k]]; nwo[k] = (nent[ord[k]] + 63) / 64; }
    for (int v = 0; v < 7; v++) lastk[v] = -1;
    for (int k = 0; k < NC; k++) { lastk[CCo[k]] = k; lastk[CDo[k]] = k; }
    for (int Xp = 0; Xp < NC; Xp++)
        for (int Y = 0; Y < NC; Y++)
            for (int s = 0; s < 2; s++)
                for (int t = 0; t < 2; t++) {
                    int r1[7] = {0}, r2[7] = {0}, v = 0;
                    r1[CCo[Xp]] = 1; r1[CDo[Xp]] = s ? -1 : 1;
                    r2[CCo[Y]] = 1; r2[CDo[Y]] = t ? -1 : 1;
                    for (int i = 0; i < 7; i++) v += r1[i] * r2[i];
                    rrt[Xp][Y][s][t] = v;
                }
    /* budget constraints: witness pairs whose contribution is not constant over the entries of some cell */
    static int *cp_tmp; static int cpt_cap = 0;
    if (cpt_cap < NPW + 1) { cpt_cap = NPW + 1; cp_tmp = realloc(cp_tmp, sizeof(int) * cpt_cap); }
    ncons = 0;
    for (int p = 0; p < NPW; p++) {
        int varies = 0, sum = 0;
        for (int x = 0; x < NC; x++) {
            int first = 1, vx = 0;
            for (int k = 0; k < ncellp[x]; k++) {
                int q = cellp[x][k].q;
                const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
                int v = va[TPa[p]] * va[TPb[p]] + vb[TPa[p]] * vb[TPb[p]];
                if (first) { vx = v; first = 0; } else if (v != vx) { varies = 1; break; }
            }
            sum += vx;
            if (varies) break;
        }
        if (varies) { if (ncons < 4096) cp_tmp[ncons++] = p; }
        else if (sum != targetAll[p]) { sv_root++; return 0; }
    }
    cpair = realloc(cpair, sizeof(int) * (ncons + 1)); ctarget = realloc(ctarget, sizeof(int) * (ncons + 1));
    for (int i = 0; i < ncons; i++) { cpair[i] = cp_tmp[i]; ctarget[i] = targetAll[cp_tmp[i]]; }
    for (int x = 0; x < NC; x++) {
        econ[x] = realloc(econ[x], (size_t)nent[x] * (ncons + 1) + 1);
        for (int t = 0; t < nent[x]; t++) {
            int q = ent[x][t].q;
            const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
            for (int i = 0; i < ncons; i++) {
                int p = cpair[i];
                econ[x][(size_t)t * ncons + i] = (signed char)(va[TPa[p]] * va[TPb[p]] + vb[TPa[p]] * vb[TPb[p]]);
            }
        }
    }
    csmin = realloc(csmin, sizeof(int) * (NC + 1) * (ncons + 1)); csmax = realloc(csmax, sizeof(int) * (NC + 1) * (ncons + 1));
    for (int i = 0; i < ncons; i++) { csmin[NC * ncons + i] = 0; csmax[NC * ncons + i] = 0; }
    for (int k = NC - 1; k >= 0; k--) {
        int x = ord[k];
        for (int i = 0; i < ncons; i++) {
            int mn = 1 << 20, mx = -(1 << 20);
            for (int t = 0; t < nent[x]; t++) { int v = econ[x][(size_t)t * ncons + i]; if (v < mn) mn = v; if (v > mx) mx = v; }
            csmin[k * ncons + i] = csmin[(k + 1) * ncons + i] + mn;
            csmax[k * ncons + i] = csmax[(k + 1) * ncons + i] + mx;
        }
    }
    for (int i = 0; i < ncons; i++) {
        crun[0][i] = 0;
        if (csmin[i] > ctarget[i] || csmax[i] < ctarget[i]) { sv_root++; return 0; }
    }
    for (int L = 0; L <= NC; L++)
        for (int Y = 0; Y < NC; Y++) {
            lev[L][Y] = realloc(lev[L][Y], sizeof(u64) * (nwo[Y] + 1));
            memset(lev[L][Y], 0, sizeof(u64) * (nwo[Y] + 1));
        }
    for (int Y = 0; Y < NC; Y++) for (int t = 0; t < nent[ord[Y]]; t++) lev[0][Y][t >> 6] |= 1ULL << (t & 63);
    /* root bounds propagation on the budgets (fixpoint) */
    if (ncons) {
        static int *bmn, *bmx; static int bcap = 0;
        if (bcap < NC * ncons) { bcap = NC * ncons; bmn = realloc(bmn, sizeof(int) * bcap); bmx = realloc(bmx, sizeof(int) * bcap); }
        for (int pass = 0; pass < 50; pass++) {
            for (int i = 0; i < ncons; i++) { dlo[i] = 0; dhi[i] = 0; }
            for (int Y = 0; Y < NC; Y++) {
                int y = ord[Y]; int *mn = bmn + Y * ncons, *mx = bmx + Y * ncons;
                for (int i = 0; i < ncons; i++) { mn[i] = 1 << 20; mx[i] = -(1 << 20); }
                int any = 0;
                for (int w = 0; w < nwo[Y]; w++) {
                    u64 xx = lev[0][Y][w];
                    while (xx) {
                        int bit = __builtin_ctzll(xx); xx &= xx - 1; any = 1;
                        const signed char *ec = econ[y] + (size_t)(64 * w + bit) * ncons;
                        for (int i = 0; i < ncons; i++) { if (ec[i] < mn[i]) mn[i] = ec[i]; if (ec[i] > mx[i]) mx[i] = ec[i]; }
                    }
                }
                if (!any) { sv_root++; return 0; }
                for (int i = 0; i < ncons; i++) { dlo[i] += mn[i]; dhi[i] += mx[i]; }
            }
            for (int i = 0; i < ncons; i++) if (dlo[i] > ctarget[i] || dhi[i] < ctarget[i]) { sv_root++; return 0; }
            int changed = 0;
            for (int Y = 0; Y < NC; Y++) {
                int y = ord[Y]; int *mn = bmn + Y * ncons, *mx = bmx + Y * ncons;
                for (int w = 0; w < nwo[Y]; w++) {
                    u64 xx = lev[0][Y][w];
                    while (xx) {
                        int bit = __builtin_ctzll(xx); xx &= xx - 1;
                        const signed char *ec = econ[y] + (size_t)(64 * w + bit) * ncons;
                        for (int i = 0; i < ncons; i++)
                            if (ec[i] + dlo[i] - mn[i] > ctarget[i] || ec[i] + dhi[i] - mx[i] < ctarget[i]) {
                                lev[0][Y][w] &= ~(1ULL << bit); changed = 1; break;
                            }
                    }
                }
            }
            if (!changed) break;
        }
        /* static suffix bounds from the propagated domains */
        for (int i = 0; i < ncons; i++) { csmin[NC * ncons + i] = 0; csmax[NC * ncons + i] = 0; }
        for (int k = NC - 1; k >= 0; k--)
            for (int i = 0; i < ncons; i++) { csmin[k * ncons + i] = csmin[(k + 1) * ncons + i] + bmn[k * ncons + i]; csmax[k * ncons + i] = csmax[(k + 1) * ncons + i] + bmx[k * ncons + i]; }
    }
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    snodes = 0;
    memset(dprof, 0, sizeof dprof);
    int r = solve_rec(0, R);
    sv_nodes_tot += snodes;
    if (getenv("DEBUGVL")) {
        fprintf(stderr, "VL ncons %d sizes", ncons);
        for (int k = 0; k < NC; k++) fprintf(stderr, " %d(%d,%d)", nent[ord[k]], CCo[k], CDo[k]);
        fprintf(stderr, " | depth");
        for (int k = 0; k <= NC; k++) fprintf(stderr, " %lld", dprof[k]);
        fputc(10, stderr);
    }
    if (r == 1) {
        sv_found++;
        if (ffound) {
            fprintf(ffound, "FOUND");
            for (int e = 0; e < 6; e++) fprintf(ffound, " %d", LG[e]);
            fputc(10, ffound);
            for (int k = 0; k < NC; k++) {
                const ENT *E = &ent[ord[k]][schosen[k]];
                const signed char *A = X + (size_t)E->va * 15, *B = X + (size_t)E->vb * 15;
                fprintf(ffound, "%d %d", CCo[k], CDo[k]);
                for (int i = 0; i < 15; i++) fprintf(ffound, " %d", E->sa * A[i]);
                for (int i = 0; i < 15; i++) fprintf(ffound, " %d", E->sb * B[i]);
                fputc(10, ffound);
            }
            fflush(ffound);
        }
    } else if (r == 2) {
        sv_limit++;
        if (flim) { for (int e = 0; e < 6; e++) fprintf(flim, "%d ", LG[e]); fputc(10, flim); fflush(flim); }
    } else sv_excl++;
    return r;
}

/* ------------------------------------------------------------------ search */
static int full_leaf(void) {
    if (!leaf_check()) return 0;
    pass_tab++;
    if (!exact_check()) return 0;
    vector_level();
    return 1;
}

static void rec(int depth) {
    nodes++;
    u64 *cand = candl[depth];
    int need = 5 - depth;
    int usestab = nroots && depth >= 1 && depth <= STABMAX;
    if (usestab) build_orbits(depth);
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
                if (usestab && ufind(uf[depth], i) != i) { stab_pruned++; continue; }
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
                full_leaf();
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

static void walk(int mode) {
    double mult = 1.0;
    for (int s = 1; s < nspan; s++) spanmark[spanel[s]] = 0;
    nspan = 1; spanel[0] = 0; spanmark[0] = 1;
    for (int depth = 0; depth < 6; depth++) {
        est_nodes += mult;
        u64 *cand = candl[depth];
        int nk = 0;
        int usestab = nroots && depth >= 1 && depth <= STABMAX;
        if (usestab) build_orbits(depth);
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
                    if (usestab && ufind(uf[depth], i) != i) continue;
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
    if (argc < 8) { puts("usage: wfe_c.exe input est|run nshards shard outprefix walks seed [pmax] [nodecap]"); return 1; }
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
    kall = malloc(nw + 1); wn = malloc(nw + 1);
    if (fread(kall, 1, nw, f) != (size_t)nw || fread(wn, 1, nw, f) != (size_t)nw) return 1;
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
    if (fread(&nroots, 4, 1, f) != 1) { puts("ERR roots"); return 1; }
    RB = malloc(sizeof(int) * (nroots + 1)); RF = malloc(sizeof(int) * (nroots + 1));
    for (int q = 0; q < nroots; q++) if (fread(&RB[q], 4, 1, f) != 1 || fread(&RF[q], 4, 1, f) != 1) { puts("ERR roots"); return 1; }
    l1off = malloc(sizeof(int) * (nvec + 1));
    long long l1cap = 1 << 20, l1n = 0;
    l1m = malloc(2 * l1cap);
    for (int v = 0; v < nvec; v++) {
        int cnt;
        l1off[v] = (int)l1n;
        if (fread(&cnt, 4, 1, f) != 1) { puts("ERR l1"); return 1; }
        while (l1n + cnt > l1cap) { l1cap *= 2; l1m = realloc(l1m, 2 * l1cap); }
        if (cnt && fread(l1m + l1n, 2, cnt, f) != (size_t)cnt) { puts("ERR l1"); return 1; }
        l1n += cnt;
    }
    l1off[nvec] = (int)l1n;
    fclose(f);
    for (int d = 0; d < 7; d++) uf[d] = malloc(sizeof(int) * (nU + 1));
    if (getenv("STABMAX")) STABMAX = atoi(getenv("STABMAX"));
    if (getenv("DYNB")) dynb = atoi(getenv("DYNB"));
    if (getenv("PAIRAC")) pair_ac = atoi(getenv("PAIRAC"));
    nshards = atoi(argv[3]); shard = atoi(argv[4]);
    int walks = atoi(argv[6]);
    rng ^= (unsigned long long)atoll(argv[7]) * 0x9E3779B97F4A7C15ULL;
    if (argc > 8) pmax = atoi(argv[8]);
    if (argc > 9) nodecap = atoll(argv[9]);
    /* L^T for the projection: GK = L L^T */
    {
        double L[15][15]; memset(L, 0, sizeof L);
        for (int i = 0; i < 15; i++)
            for (int j = 0; j <= i; j++) {
                double s = GK[i][j];
                for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
                if (i == j) { if (s <= 0) { puts("ERR chol"); return 1; } L[i][i] = sqrt(s); }
                else L[i][j] = s / L[j][j];
            }
        for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) Lt[i][j] = L[j][i];
    }
    GWV = malloc(sizeof(int) * 15 * (size_t)nvec); LWV = malloc(sizeof(double) * 15 * (size_t)nvec);
    for (int v = 0; v < nvec; v++)
        for (int i = 0; i < 15; i++) {
            int s1 = 0; double t1 = 0;
            for (int j = 0; j < 15; j++) { s1 += GK[i][j] * X[(size_t)v * 15 + j]; t1 += Lt[i][j] * X[(size_t)v * 15 + j]; }
            GWV[(size_t)v * 15 + i] = s1; LWV[(size_t)v * 15 + i] = t1;
        }
    /* value masks per vector */
    M26 = calloc(nvec + 1, 8); M37 = calloc(nvec + 1, 8);
    for (int v = 0; v < nvec; v++)
        for (int t = 0; t < nw; t++) {
            int a = VAL[(size_t)v * nw + t]; if (a < 0) a = -a;
            if (wn[t] == 6 && a == 2) M26[v] |= 1ULL << t;
            if (wn[t] == 7 && a == 3) M37[v] |= 1ULL << t;
        }
    vstart = malloc(sizeof(int) * nU); vend = malloc(sizeof(int) * nU);
    for (int i = 0; i < nU; i++) { vstart[i] = 0; vend[i] = 0; }
    for (int v = 0; v < nvec; v++) { int c = vcls[v]; if (v == 0 || vcls[v - 1] != c) vstart[c] = v; vend[c] = v + 1; }
    int *GW = malloc(sizeof(int) * 15 * (size_t)nvec);
    for (int v = 0; v < nvec; v++)
        for (int k = 0; k < 15; k++) { int s = 0; for (int j = 0; j < 15; j++) s += GK[k][j] * X[(size_t)v * 15 + j]; GW[(size_t)v * 15 + k] = s; }
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
    free(GW);
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
        for (long long q = pst[c]; q < pen[c]; q++)
            for (int s = 0; s < 2; s++) {
                int16_t m = s ? PM2[q] : PM1[q];
                int a = mst[c], b = men[c] - 1;
                while (a < b) { int mid = (a + b) / 2; if (UM[mid] < m) a = mid + 1; else b = mid; }
                if (s) R2[q] = a; else R1[q] = a;
            }
    /* pairs indexed by a-mask row */
    RPoff = calloc((size_t)nrows + 2, sizeof(long long));
    for (long long q = 0; q < npairs; q++) { RPoff[R1[q] + 1]++; RPoff[R2[q] + 1]++; }
    for (long long r = 0; r < nrows; r++) RPoff[r + 1] += RPoff[r];
    RPq = malloc(sizeof(int) * (size_t)(2 * npairs + 1)); RPo = malloc((size_t)(2 * npairs + 1));
    {
        long long *fill = malloc(sizeof(long long) * (size_t)(nrows + 1));
        memcpy(fill, RPoff, sizeof(long long) * (size_t)(nrows + 1));
        for (long long q = 0; q < npairs; q++) {
            RPq[fill[R1[q]]] = (int)q; RPo[fill[R1[q]]++] = 1;
            RPq[fill[R2[q]]] = (int)q; RPo[fill[R2[q]]++] = 2;
        }
        free(fill);
    }
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
    fprintf(stderr, "vectors %d classes %d pairs %lld mask rows %lld witnesses %d (pairs %d) roots %d l1 %lld\n", nvec, nU,
            npairs, nrows, nw, NPW, nroots, l1n);
    killcnt = calloc(NPW + 1, sizeof(long long));
    alo = malloc(sizeof(int) * (NPW + 1)); ahi = malloc(sizeof(int) * (NPW + 1));
    acmn = malloc(sizeof(int) * (NPW + 1)); acmx = malloc(sizeof(int) * (NPW + 1));
    int maxs = walks > 0 ? walks : 1;
    unsigned char *ks = calloc((size_t)maxs * NPW, 1);
    long long *prev = calloc(NPW + 1, sizeof(long long));
    int nks = 0;
    for (int k = 0; k < walks; k++) {
        long long before = samp_killed;
        memcpy(prev, killcnt, sizeof(long long) * NPW);
        walk(0);
        if (samp_killed > before) {
            for (int p = 0; p < NPW; p++) if (killcnt[p] > prev[p]) ks[(size_t)nks * NPW + p] = 1;
            nks++;
        }
    }
    fprintf(stderr, "sampling: walks %d sampled leaves %lld killed %lld m1-dead %lld | est leaves %.4g est M1-valid %.4g\n",
            walks, samp_leaves, samp_killed, samp_m1, est_leaves / (walks ? walks : 1), est_m1 / (walks ? walks : 1));
    Plist = malloc(sizeof(int) * (NPW + 1)); nP = 0;
    unsigned char *inP = calloc(NPW + 1, 1);
    for (int rnd = 0; rnd < 2; rnd++)
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
    TMN = malloc((size_t)nrows * nP + 1); TMX = malloc((size_t)nrows * nP + 1);
    for (long long r = 0; r < nrows * nP; r++) { TMN[r] = 127; TMX[r] = -128; }
    for (int c = 0; c < nU; c++)
        for (long long q = pst[c]; q < pen[c]; q++) {
            const signed char *va = VAL + (size_t)PI[q] * nw, *vb = VAL + (size_t)PJ[q] * nw;
            for (int p = 0; p < nP; p++) {
                int P_ = Plist[p];
                int v = va[TPa[P_]] * va[TPb[P_]] + vb[TPa[P_]] * vb[TPb[P_]];
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
        for (int k = 0; k < walks; k++) walk(1);
        printf("KNUTH walks %d est_nodes %.4g est_leaves %.4g est_survivors %.4g\n", walks, est_nodes / walks,
               est_leaves / walks, est_surv / walks);
        return 0;
    }
    if (!strcmp(argv[2], "check")) {          /* check mode: skeletons from argv[10], one per line (6 class values) */
        FILE *g = fopen(argv[10], "r");
        int v[6];
        while (g && fscanf(g, "%d %d %d %d %d %d", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]) == 6) {
            int ok = 1;
            for (int e = 0; e < 6; e++) { chosen[e] = idxof[v[e]]; if (chosen[e] < 0) ok = 0; }
            if (!ok) { puts("CHECK notusable"); continue; }
            long long e0 = ex_empty, b0 = ex_budget, r0 = sv_root, x0 = sv_excl, l0 = sv_limit, f0 = sv_found, n0 = sv_nodes_tot;
            long long m0 = m1dead, t0 = bdead, a0 = ex_ac;
            int tab = leaf_check();
            int ex = exact_check();
            const char *vl = "-";
            if (ex) {
                vector_level();
                vl = sv_root > r0 ? "root" : sv_excl > x0 ? "excluded" : sv_limit > l0 ? "LIMIT" : sv_found > f0 ? "FOUND" : "?";
            }
            printf("CHECK tab=%d(%s) exact=%s vl=%s nodes=%lld%c", tab, m1dead > m0 ? "m1" : bdead > t0 ? "budget" : "pass",
                   ex ? "pass" : (ex_empty > e0 ? "empty" : (ex_ac > a0 ? "ac" : "budget")), vl, sv_nodes_tot - n0, 10);
        }
        return 0;
    }
    char fn[1024];
    snprintf(fn, sizeof fn, "%s.limit%d", argv[5], shard); flim = fopen(fn, "w");
    snprintf(fn, sizeof fn, "%s.found%d", argv[5], shard); ffound = fopen(fn, "w");
    memset(spanmark, 0, sizeof spanmark); spanmark[0] = 1; spanel[0] = 0; nspan = 1;
    memset(ocnt, 0, sizeof ocnt); memset(ucnt, 0, sizeof ucnt);
    rec(0);
    printf("DONE nodes %lld leaves %lld tab_m1 %lld tab_budget %lld tab_pass %lld ex_empty %lld ex_budget %lld sv_root %lld "
           "sv_excluded %lld sv_limit %lld sv_found %lld sv_nodes %lld stab_pruned %lld ex_ac %lld\n", nodes, leaves, m1dead, bdead,
           pass_tab, ex_empty, ex_budget, sv_root, sv_excl, sv_limit, sv_found, sv_nodes_tot, stab_pruned, ex_ac);
    fclose(flim); fclose(ffound);
    return 0;
}
