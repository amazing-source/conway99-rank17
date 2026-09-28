/* Full V = Z^17 pipeline in C (structures, 28/09): stage A (as stageA.c: partition refinement, S_7 breaking, usability,
   incremental M1, PE, coverage) + per-leaf stage B (domains as stageB.Skeleton.build: refined L1, M1 keyed on (w.pi_e) mod 4,
   sign normalisation; then exact DFS: plain with a small cap, then bitsets + 2-open star look-ahead as sbdfs3.c).
   Input: the text file of mk_stageA_input.py.
   Usage: z17pipe.exe INPUT OUTPREFIX [sample_prob seed] [part r m] [heavycap]
   Output: OUTPREFIX.stats (final JSON-like line) and OUTPREFIX.surv (FOUND / LIMIT leaves with columns).                   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
#define NC 21

/* ---------------- globals: pair data ---------------- */
static int Y1[17], Y2[17], O1 = 0, O2 = 0;
static int npe, peB[4096], peK[4096], peS[4096];
static unsigned char *usab;
static int **amask, *namask;
static signed char **cvec; static int *ncvec;
static double Bw[17 * 15];                 /* orthonormal basis of W (row-major 17x15) */
static FILE *fsurv;
static double samp = 1.0; static u64 rng_state = 88172645463325252ULL;
static int part_r = 0, part_m = 1; static long long top_counter = 0;
static long long heavycap = 200000000LL;
static long long st_leaves = 0, st_proc = 0, st_empty = 0, st_root = 0, st_deeper = 0, st_heavy = 0, st_found = 0,
                 st_limit = 0, m1_prunes = 0, st_nodes = 0;
static time_t t_start, t_last;
static long long debugN = 0; static FILE *fdbg = NULL;
static double tm_build = 0, tm_prep = 0, tm_solve = 0;
static int killer[NC + 1];
static int countonly = 0, toponly = 0; static long long top_total = 0;
static double probe_q = 1.0;
static long maxtime = 0; static int stop_now = 0;

static int pc(unsigned x) { return __builtin_popcount(x); }
static double urand(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return (rng_state >> 11) * (1.0 / 9007199254740992.0); }
static int canon7(int x) { return (x & 1) ? (x ^ 127) : x; }

/* ---------------- class vectors, usability, abar masks ---------------- */
static int usable(unsigned S) {
    if (usab[S]) return usab[S] == 1;
    int idx[17], n = 0, rest[17], nr = 0;
    for (int i = 0; i < 17; i++) { if (S >> i & 1) idx[n++] = i; else rest[nr++] = i; }
    static signed char V[20000][17];
    int nv = 0, ntw = (n == 10) ? 0 : (n == 6 ? 1 : (n == 2 ? 2 : -1));
    if (ntw < 0) { usab[S] = 2; return 0; }
    for (int sg = 0; sg < (1 << n); sg++) {
        signed char w[17];
        memset(w, 0, 17);
        for (int a = 0; a < n; a++) w[idx[a]] = (sg >> a & 1) ? -1 : 1;
#define PUSHIF { int d1 = 0, d2 = 0; for (int i = 0; i < 17; i++) { d1 += w[i] * Y1[i]; d2 += w[i] * Y2[i]; } if (!d1 && !d2) memcpy(V[nv++], w, 17); }
        if (ntw == 0) PUSHIF
        else if (ntw == 1) {
            for (int a = 0; a < nr; a++) for (int s = -2; s <= 2; s += 4) { w[rest[a]] = (signed char)s; PUSHIF w[rest[a]] = 0; }
        } else {
            for (int a = 0; a < nr; a++) for (int b = a + 1; b < nr; b++)
                for (int s = -2; s <= 2; s += 4) for (int t = -2; t <= 2; t += 4) {
                    w[rest[a]] = (signed char)s; w[rest[b]] = (signed char)t; PUSHIF w[rest[a]] = 0; w[rest[b]] = 0;
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
    amask[S] = malloc(sizeof(int) * (size_t)(nl ? nl : 1)); memcpy(amask[S], list, sizeof(int) * (size_t)nl); namask[S] = nl;
    cvec[S] = malloc((size_t)(nv ? nv : 1) * 17); memcpy(cvec[S], V, (size_t)nv * 17); ncvec[S] = nv;
    usab[S] = ok ? 1 : 2;
    return ok;
}

/* ---------------- stage B: domains and DFS ---------------- */
static int CC[NC], CD[NC], ndom[NC], nwd[NC];
static const signed char *cv[NC]; static int ncv[NC];     /* class vectors of the cell (shared cache) */
static int *ipair[NC];                                    /* 2 indices per domain entry into cv[X] */
static int dotbuf[2][20000];
static signed char *dom[NC];
static double *proj[NC];
static u64 *kp[NC], *km[NC];
static int chosen[NC], rrt[NC][NC][2][2], lastk[7];
static long long nodes, cap;
static int *levi[NC + 1][NC], levn[NC + 1][NC];     /* plain solver: index lists */
static u64 *levb[NC + 1][NC];                       /* bitset solver */
static u64 *compat[NC][NC];
static unsigned char *rowdone[NC];
static void compat_row(int X, int t) {
    const signed char *a = dom[X] + 34 * t;
    for (int Y = X + 1; Y < NC; Y++) {
        u64 *row = compat[X][Y] + (size_t)t * nwd[Y];
        for (int q = 0; q < ndom[Y]; q++) {
            const signed char *b = dom[Y] + 34 * q;
            int ok = 1;
            for (int s2 = 0; s2 < 2 && ok; s2++) for (int tt = 0; tt < 2 && ok; tt++) {
                int g = rrt[X][Y][s2][tt];
                for (int i = 0; i < 17; i++) g += a[17 * s2 + i] * b[17 * tt + i];
                if (!(g == 0 || g == 4 || g == -4)) ok = 0;
            }
            if (ok) row[q >> 6] |= 1ULL << (q & 63);
        }
    }
    rowdone[X][t] = 1;
}
static u64 *htab; static int hcap = 0;
/* 3-open star look-ahead: per vertex v, its last three cells L[v][0..2]; pair sum-sets H[v][pair] of contrib_v(p)+contrib_v(q)
   over M2-compatible (p,q) of the ORIGINAL domains (a superset of what is still possible: valid necessary condition).
   Keys: 17 entries in [-8,8] packed 5 bits each into two u64 (lo: coords 0..11, hi: 12..16). */
static int L3v[7][3], has3[7][3];
static u64 *H3[7][3]; static long long H3n[7][3];
static int anybit(const u64 *b, int n);
static long long st_sum3 = 0, st_acwipe = 0;
static int use_ac = 1;
/* compat[Y][X] for Y > X is not stored; revision of X by Y uses compat[X][Y] (X<Y) or the transpose test */
static u64 *cmT[NC][NC];
static int ac_bits(u64 **al) {
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int X = 0; X < NC; X++)
            for (int Y = 0; Y < NC; Y++) {
                if (X == Y) continue;
                for (int w = 0; w < nwd[X]; w++) {
                    u64 x = al[X][w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        int t = 64 * w + bit, sup = 0;
                        if (X < Y) {
                            const u64 *row = compat[X][Y] + (size_t)t * nwd[Y];
                            for (int i = 0; i < nwd[Y]; i++) if (row[i] & al[Y][i]) { sup = 1; break; }
                        } else {
                            const u64 *row = cmT[X][Y] + (size_t)t * nwd[Y];
                            for (int i = 0; i < nwd[Y]; i++) if (row[i] & al[Y][i]) { sup = 1; break; }
                        }
                        if (!sup) { al[X][w] &= ~(1ULL << bit); changed = 1; }
                    }
                }
                if (!anybit(al[X], nwd[X])) return 0;
            }
    }
    return 1;
}
static long long st_zsac_codes = 0, st_zsac_acdead = 0, st_zsac_search = 0;
static int use_zsac = 1;
static int ac_full(void) {
    for (int X = 0; X < NC; X++) for (int t = 0; t < ndom[X]; t++) if (!rowdone[X][t]) compat_row(X, t);
    u64 *al[NC];
    for (int X = 0; X < NC; X++) al[X] = levb[0][X];
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int X = 0; X < NC; X++)
            for (int Y = 0; Y < NC; Y++) {
                if (X == Y) continue;
                for (int w = 0; w < nwd[X]; w++) {
                    u64 x = al[X][w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        int t = 64 * w + bit, sup = 0;
                        if (X < Y) {
                            const u64 *row = compat[X][Y] + (size_t)t * nwd[Y];
                            for (int i = 0; i < nwd[Y]; i++) if (row[i] & al[Y][i]) { sup = 1; break; }
                        } else {
                            /* Y < X: t in D_X is supported iff some alive q in D_Y has bit t in compat[Y][X] row q */
                            for (int wq = 0; wq < nwd[Y] && !sup; wq++) {
                                u64 y = al[Y][wq];
                                while (y) {
                                    int bq = __builtin_ctzll(y); y &= y - 1;
                                    const u64 *row = compat[Y][X] + (size_t)(64 * wq + bq) * nwd[X];
                                    if (row[t >> 6] >> (t & 63) & 1) { sup = 1; break; }
                                }
                            }
                        }
                        if (!sup) { al[X][w] &= ~(1ULL << bit); changed = 1; }
                    }
                }
                if (!anybit(al[X], nwd[X])) return 0;
            }
    }
    return 1;
}
static long long look3_cap = 150000;
static void key85(const int *v, u64 *lo, u64 *hi) {
    u64 a = 0, b = 0;
    for (int i = 0; i < 12; i++) a |= (u64)(v[i] + 16) << (5 * i);
    for (int i = 12; i < 17; i++) b |= (u64)(v[i] + 16) << (5 * (i - 12));
    *lo = a; *hi = b | (1ULL << 40);
}
static u64 hmix(u64 a, u64 b) { u64 h = a * 0x9E3779B97F4A7C15ULL ^ b * 0xC2B2AE3D27D4EB4FULL; return h ^ (h >> 31); }
static void h3_insert(u64 *T, long long n, u64 a, u64 b) {
    long long q = (long long)(hmix(a, b) % (u64)n);
    while (T[2 * q + 1]) { if (T[2 * q] == a && T[2 * q + 1] == b) return; q = (q + 1) % n; }
    T[2 * q] = a; T[2 * q + 1] = b;
}
static int h3_has(const u64 *T, long long n, u64 a, u64 b) {
    long long q = (long long)(hmix(a, b) % (u64)n);
    while (T[2 * q + 1]) { if (T[2 * q] == a && T[2 * q + 1] == b) return 1; q = (q + 1) % n; }
    return 0;
}
static void contrib(int X, int t, int v, int *out) {
    const signed char *u = dom[X] + 34 * t;
    for (int i = 0; i < 17; i++) out[i] = (CC[X] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
}

static unsigned long long ckey(const signed char *u, int plus, const int *v) {
    unsigned long long k = 0;
    for (int i = 16; i >= 0; i--) { int c = u ? (plus ? u[i] + u[17 + i] : u[i] - u[17 + i]) : v[i]; k = k * 9ULL + (unsigned long long)(c + 4); }
    return k;
}
static int psd_ok(double R[15][15]) {
    double L[15][15];
    for (int i = 0; i < 15; i++)
        for (int j = 0; j <= i; j++) {
            double s = R[i][j] + (i == j ? 1e-6 : 0.0);
            for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
            if (i == j) { if (s <= 0) return 0; L[i][i] = sqrt(s); } else L[i][j] = s / L[j][j];
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
static int first_pos(const int *a) { for (int i = 0; i < 17; i++) if (a[i]) return a[i] > 0; return 0; }
/* symmetry C: sign flips of coordinates where y1 = y2 = 0 (zc_mask).  H-normalisation must be C-invariant, so it looks only at
   coordinates outside zc_mask; a vector with no support there is not normalised (both signs kept).  In the search, a lex-leader
   rule fixes the C-signs: along the path, the first non-zero entry at each zc coordinate (w+ before w-) must be positive. */
static unsigned zc_mask = 0;
static int use_csign = 1;
static unsigned cfix[64];
static int first_pos_nz(const int *a) {           /* 1: positive, 0: negative, 2: no support outside zc (accept) */
    for (int i = 0; i < 17; i++) if (!(zc_mask >> i & 1) && a[i]) return a[i] > 0;
    return 2;
}
static int csign_ok(const signed char *w, unsigned fixed, unsigned *out) {
    unsigned f = fixed;
    if (use_csign) {
        for (int h = 0; h < 2; h++) {
            const signed char *v = w + 17 * h;
            for (int i = 0; i < 17; i++) {
                if (!(zc_mask >> i & 1) || (f >> i & 1) || !v[i]) continue;
                if (v[i] < 0) return 0;
                f |= 1u << i;
            }
        }
    }
    *out = f;
    return 1;
}

static const int *g_keys;
static int cmpk(const void *a, const void *b) { int x = g_keys[*(const int *)a], y = g_keys[*(const int *)b]; return (x > y) - (x < y); }

/* build domains for the current leaf; returns 0 if some cell is empty */
static int build(const unsigned *rows, const int *cols) {
    for (int X = 0; X < NC; X++) {
        int c = CC[X], d = CD[X];
        unsigned S = rows[c] ^ rows[d];
        int n = ncvec[S];
        const signed char *W = cvec[S];
        int odd[17], no = 0;
        for (int i = 0; i < 17; i++) if (S >> i & 1) odd[no++] = i;
        static unsigned char keepv[20000];
        for (int t = 0; t < n; t++) keepv[t] = 1;
        if (no == 2) {
            int i = odd[0], j = odd[1];
            int forb[4] = {0, canon7(1 << c), canon7(1 << d), canon7((1 << c) ^ (1 << d))};
            int span[4] = {0, cols[i], cols[j], cols[i] ^ cols[j]};
            for (int t = 0; t < n; t++) {
                const signed char *w = W + 17 * t;
                int kl[2], nk = 0;
                for (int x = 0; x < 17; x++) if (w[x] == 2 || w[x] == -2) kl[nk++] = x;
                int base = cols[kl[0]] ^ cols[kl[1]], bad = 0;
                for (int s = 0; s < 4 && !bad; s++) { int v = canon7(base ^ span[s]); for (int f = 0; f < 4; f++) if (v == forb[f]) bad = 1; }
                if (bad) keepv[t] = 0;
            }
        }
        /* keys */
        static int key[20000], need[20000], ord[20000];
        for (int t = 0; t < n; t++) {
            const signed char *w = W + 17 * t;
            int kk = 0, nn = 0, p4 = 1;
            for (int e = 0; e < 7; e++) {
                int s = 0;
                for (int i = 0; i < 17; i++) if (rows[e] >> i & 1) s += w[i];
                int te = ((s % 4) + 4) % 4;
                int tg = (c == 0) ? (e != 0) : (e == c);
                kk += te * p4; nn += (((2 * tg - te) % 4) + 4) % 4 * p4; p4 *= 4;
            }
            key[t] = kk; need[t] = nn;
        }
        for (int t = 0; t < n; t++) ord[t] = t;
        g_keys = key;
        qsort(ord, (size_t)n, sizeof(int), cmpk);
        /* collect pairs */
        int capd = 1024, nd = 0;
        signed char *D = malloc((size_t)capd * 34);
        int *IP = malloc(sizeof(int) * 2 * (size_t)capd);
        for (int p = 0; p < n; p++) {
            if (!keepv[p]) continue;
            int nk = need[p];
            /* binary search range of key == nk in ord */
            int lo = 0, hi = n;
            while (lo < hi) { int mid = (lo + hi) / 2; if (key[ord[mid]] < nk) lo = mid + 1; else hi = mid; }
            for (int r = lo; r < n && key[ord[r]] == nk; r++) {
                int q = ord[r];
                if (q == p || !keepv[q]) continue;
                const signed char *wp = W + 17 * p, *wq = W + 17 * q;
                int g = 0;
                for (int i = 0; i < 17; i++) g += wp[i] * wq[i];
                if (!(g == 0 || g == 4 || g == -4)) continue;
                if (c == 0) {
                    int a[17], b[17];
                    for (int i = 0; i < 17; i++) { a[i] = (wp[i] + wq[i]) / 2; b[i] = (wp[i] - wq[i]) / 2; }
                    if (use_csign) {
                        if (first_pos_nz(b) == 0) continue;
                        if (d == 1 && first_pos_nz(a) == 0) continue;
                    } else {
                        if (!first_pos(b)) continue;
                        if (d == 1 && !first_pos(a)) continue;
                    }
                }
                if (nd == capd) { capd *= 2; D = realloc(D, (size_t)capd * 34); IP = realloc(IP, sizeof(int) * 2 * (size_t)capd); }
                memcpy(D + 34 * nd, wp, 17); memcpy(D + 34 * nd + 17, wq, 17);
                IP[2 * nd] = p; IP[2 * nd + 1] = q; nd++;
            }
        }
        dom[X] = D; ndom[X] = nd; ipair[X] = IP; cv[X] = W; ncv[X] = n;
        if (!nd) {
            for (int Z = 0; Z <= X; Z++) { free(dom[Z]); dom[Z] = NULL; free(ipair[Z]); ipair[Z] = NULL; }
            return 0;
        }
    }
    return 1;
}

/* exact mod-4 lift test (see z4test.py): self-dual Z4 codes with residue D <-> 6x6 matrix M over F2 (21 free bits) */
static int use_z4 = 1;
static long long st_z4 = 0;
static unsigned bz[7];          /* bz[c], c=1..6: 24-bit basis masks (bits 0..6 e-part, 7..23 f-part) */
static int Jz[7][7], kapz[7][7];
static u64 SigZ[NC];            /* per cell: 64-bit set of 6-bit signatures */
static unsigned char *sigE[NC];  /* per entry signature */
static int use_zfilter = 1;
static long long st_zfilt_in = 0, st_zfilt_out = 0, st_zM = 0;
static u64 allowZ[NC];
static int Mz[7];
static int cellidx[7][7];
static int z4_rec(int c) {
    if (c == 7) return 1;
    u64 S = SigZ[cellidx[0][c]];
    while (S) {
        int m = __builtin_ctzll(S); S &= S - 1;
        int ok = 1;
        for (int e = 1; e < c && ok; e++) if ((((m >> (e - 1)) & 1) ^ ((Mz[e] >> (c - 1)) & 1)) != Jz[c][e]) ok = 0;
        for (int e = 1; e < c && ok; e++) { int v = Mz[e] ^ m ^ kapz[e][c]; if (!(SigZ[cellidx[e][c]] >> v & 1)) ok = 0; }
        if (!ok) continue;
        Mz[c] = m;
        if (z4_rec(c + 1)) return 1;
    }
    return 0;
}
static int z4_ok(const unsigned *rows) {
    for (int c = 1; c <= 6; c++) bz[c] = 1u | (1u << c) | (rows[c] << 7);
    for (int c = 1; c <= 6; c++) for (int e = 1; e <= 6; e++) Jz[c][e] = (pc(bz[c] & bz[e]) / 2) & 1;
    for (int c = 1; c <= 6; c++) for (int d = 1; d <= 6; d++) {
        int v = 0;
        for (int e = 1; e <= 6; e++) if (pc(bz[c] & bz[d] & bz[e]) & 1) v |= 1 << (e - 1);
        kapz[c][d] = v;
    }
    for (int X = 0; X < NC; X++) {
        cellidx[CC[X]][CD[X]] = X;
        free(sigE[X]); sigE[X] = malloc((size_t)(ndom[X] ? ndom[X] : 1));
        u64 S = 0;
        int c = CC[X], d = CD[X];
        for (int t = 0; t < ndom[X]; t++) {
            const signed char *w = dom[X] + 34 * t;           /* w+ ; q+ = (e_c + e_d | w+) */
            unsigned x = 0;                                   /* e-part: +1 entries give x = 0 */
            for (int i = 0; i < 17; i++) {
                int q = w[i], xi;
                if (q & 1) xi = ((((q - 1) % 4) + 4) % 4) / 2; else xi = (((q % 4) + 4) % 4) / 2;
                if (xi) x |= 1u << (7 + i);
            }
            int sig = 0;
            for (int e = 1; e <= 6; e++) if (pc(x & bz[e]) & 1) sig |= 1 << (e - 1);
            S |= 1ULL << sig;
            sigE[X][t] = (unsigned char)sig;
        }
        SigZ[X] = S;
        (void)c; (void)d;
    }
    return z4_rec(1);
}

/* Z4 code-consistency filter: allowed signatures per cell = targets over ALL feasible M (sound: every solution lies in one Q_M) */
static int zfM[7];
#define ZMCAP 16384
static unsigned char *zT = NULL;      /* nZM x NC targets */
static int nZM = 0, zmOverflow = 0;
static void zf_rec(int c) {
    if (c == 7) {
        st_zM++;
        if (nZM < ZMCAP) {
            if (!zT) zT = malloc((size_t)ZMCAP * NC);
            for (int X = 0; X < NC; X++) {
                int cc = CC[X], dd = CD[X];
                zT[(size_t)nZM * NC + X] = (unsigned char)((cc == 0) ? zfM[dd] : (zfM[cc] ^ zfM[dd] ^ kapz[cc][dd]));
            }
            nZM++;
        } else zmOverflow = 1;
        for (int X = 0; X < NC; X++) {
            int cc = CC[X], dd = CD[X];
            int tgt = (cc == 0) ? zfM[dd] : (zfM[cc] ^ zfM[dd] ^ kapz[cc][dd]);
            allowZ[X] |= 1ULL << tgt;
        }
        return;
    }
    u64 S = SigZ[cellidx[0][c]];
    while (S) {
        int m = __builtin_ctzll(S); S &= S - 1;
        int ok = 1;
        for (int e = 1; e < c && ok; e++) if ((((m >> (e - 1)) & 1) ^ ((zfM[e] >> (c - 1)) & 1)) != Jz[c][e]) ok = 0;
        for (int e = 1; e < c && ok; e++) { int v = zfM[e] ^ m ^ kapz[e][c]; if (!(SigZ[cellidx[e][c]] >> v & 1)) ok = 0; }
        if (!ok) continue;
        zfM[c] = m;
        zf_rec(c + 1);
    }
}
static void zfilter(void) {
    for (int X = 0; X < NC; X++) allowZ[X] = 0;
    nZM = 0; zmOverflow = 0;
    zf_rec(1);
    for (int X = 0; X < NC; X++) {
        int m = 0;
        st_zfilt_in += ndom[X];
        for (int t = 0; t < ndom[X]; t++) {
            if (!(allowZ[X] >> sigE[X][t] & 1)) continue;
            if (m != t) {
                memcpy(dom[X] + 34 * m, dom[X] + 34 * t, 34);
                ipair[X][2 * m] = ipair[X][2 * t]; ipair[X][2 * m + 1] = ipair[X][2 * t + 1];
                sigE[X][m] = sigE[X][t];
            }
            m++;
        }
        ndom[X] = m;
        st_zfilt_out += m;
    }
}

/* plain DFS (sbdfs.c) */
static int rec_plain(int k, double R[15][15]);
static int try_plain(int k, int t, double R[15][15]) {
    const signed char *w = dom[k] + 34 * t;
    if (!csign_ok(w, cfix[k], &cfix[k + 1])) return 0;
    const double *pv = proj[k] + 30 * t;
    double R2[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pv[i] * pv[j] - pv[15 + i] * pv[15 + j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    chosen[k] = t;
    int s[17];
    for (int v = 0; v < 7; v++) if (lastk[v] == k) { starsum(v, k, s); for (int i = 0; i < 17; i++) if (s[i]) return 0; }
    for (int yy = k; yy < NC; yy++) {
        int Y;
        if (yy == k) { Y = killer[k]; if (Y <= k || Y >= NC) continue; }
        else { Y = yy; if (Y == k || (Y == killer[k] && killer[k] > k)) continue; }
        int m = 0, *src = levi[k][Y], ns = levn[k][Y], *dst = levi[k + 1][Y];
        int usevec = (4 * ns > 2 * ncv[Y]);
        if (usevec) {
            for (int j = 0; j < ncv[Y]; j++) {
                const signed char *u = cv[Y] + 17 * j;
                int g0 = 0, g1 = 0;
                for (int i = 0; i < 17; i++) { g0 += w[i] * u[i]; g1 += w[17 + i] * u[i]; }
                dotbuf[0][j] = g0; dotbuf[1][j] = g1;
            }
        }
        for (int q = 0; q < ns; q++) {
            int ok = 1;
            if (usevec) {
                int ip = ipair[Y][2 * src[q]], iq = ipair[Y][2 * src[q] + 1];
                int g;
                g = dotbuf[0][ip] + rrt[k][Y][0][0]; if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                if (ok) { g = dotbuf[0][iq] + rrt[k][Y][0][1]; if (!(g == 0 || g == 4 || g == -4)) ok = 0; }
                if (ok) { g = dotbuf[1][ip] + rrt[k][Y][1][0]; if (!(g == 0 || g == 4 || g == -4)) ok = 0; }
                if (ok) { g = dotbuf[1][iq] + rrt[k][Y][1][1]; if (!(g == 0 || g == 4 || g == -4)) ok = 0; }
            } else {
                const signed char *u = dom[Y] + 34 * src[q];
                for (int a = 0; a < 2 && ok; a++) for (int b = 0; b < 2 && ok; b++) {
                    int g = rrt[k][Y][a][b];
                    for (int i = 0; i < 17; i++) g += w[17 * a + i] * u[17 * b + i];
                    if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                }
            }
            if (ok) dst[m++] = src[q];
        }
        if (!m) { killer[k] = Y; return 0; }
        levn[k + 1][Y] = m;
    }
    for (int v = 0; v < 7; v++) {
        int open = -1, cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CC[Y] == v || CD[Y] == v) { open = Y; cnt++; }
        if (cnt != 1) continue;
        starsum(v, k, s);
        int *arr = levi[k + 1][open], n = levn[k + 1][open], m = 0;
        for (int q = 0; q < n; q++) {
            const signed char *u = dom[open] + 34 * arr[q];
            int ok = 1;
            for (int i = 0; i < 17 && ok; i++) { int c = (CC[open] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i]; if (c != -s[i]) ok = 0; }
            if (ok) arr[m++] = arr[q];
        }
        if (!m) return 0;
        levn[k + 1][open] = m;
    }
    return rec_plain(k + 1, R2);
}
static int rec_plain(int k, double R[15][15]) {
    nodes++;
    if (nodes > cap) return 2;
    if (k == NC) { for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0; return 1; }
    for (int q = 0; q < levn[k][k]; q++) { int r = try_plain(k, levi[k][k][q], R); if (r) return r; }
    return 0;
}

/* bitset DFS (sbdfs3.c) */
static int anybit(const u64 *b, int n) { for (int i = 0; i < n; i++) if (b[i]) return 1; return 0; }
static int rec_bits(int k, double R[15][15]);
static int use_dom = 1;
/* cmT declared earlier: cmT[Y][X] = rows over D_X for each value of Y, Y > X (transpose of compat[X][Y]) */
static u64 *dlev[NC + 1][NC];            /* alive bitsets per depth */
static int assigned[NC];
static int rec_dom(int depth, double R[15][15]);
static const u64 *crow(int X, int Y, int t) {     /* bitset over D_Y of values compatible with X = t */
    if (X < Y) return compat[X][Y] + (size_t)t * nwd[Y];
    return cmT[X][Y] + (size_t)t * nwd[Y];
}
static int popc(const u64 *b, int n) { int c = 0; for (int i = 0; i < n; i++) c += __builtin_popcountll(b[i]); return c; }
static void starsum_dom(int v, int *s) {
    memset(s, 0, 17 * sizeof(int));
    for (int x = 0; x < NC; x++) {
        if (!assigned[x]) continue;
        const signed char *w = dom[x] + 34 * chosen[x];
        if (CC[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] + w[17 + i];
        else if (CD[x] == v) for (int i = 0; i < 17; i++) s[i] += w[i] - w[17 + i];
    }
}
static int try_dom(int depth, int X, int t, double R[15][15]) {
    if (!csign_ok(dom[X] + 34 * t, cfix[depth], &cfix[depth + 1])) return 0;
    const double *pv = proj[X] + 30 * t;
    double R2[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pv[i] * pv[j] - pv[15 + i] * pv[15 + j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    chosen[X] = t; assigned[X] = 1;
    int s[17], res = 0;
    for (int v = 0; v < 7; v++) {                        /* star closure for vertices now complete */
        if (CC[X] != v && CD[X] != v) continue;
        int open = 0;
        for (int Y = 0; Y < NC; Y++) if ((CC[Y] == v || CD[Y] == v) && !assigned[Y]) open++;
        if (!open) { starsum_dom(v, s); for (int i = 0; i < 17; i++) if (s[i]) goto out; }
    }
    for (int Y = 0; Y < NC; Y++) {
        if (assigned[Y]) continue;
        const u64 *cm = crow(X, Y, t);
        u64 *src = dlev[depth][Y], *dst = dlev[depth + 1][Y], acc = 0;
        for (int i = 0; i < nwd[Y]; i++) { dst[i] = src[i] & cm[i]; acc |= dst[i]; }
        if (!acc) goto out;
    }
    for (int v = 0; v < 7; v++) {                        /* star look-ahead with 1 or 2 open cells */
        int op[2], cnt = 0;
        for (int Y = 0; Y < NC; Y++) if ((CC[Y] == v || CD[Y] == v) && !assigned[Y]) { if (cnt < 2) op[cnt] = Y; cnt++; }
        if (cnt == 1) {
            int Y = op[0];
            starsum_dom(v, s);
            u64 *b = dlev[depth + 1][Y];
            for (int w = 0; w < nwd[Y]; w++) {
                u64 x = b[w];
                while (x) {
                    int bit = __builtin_ctzll(x); x &= x - 1;
                    const signed char *u = dom[Y] + 34 * (64 * w + bit);
                    for (int i = 0; i < 17; i++) { int c = (CC[Y] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i]; if (c != -s[i]) { b[w] &= ~(1ULL << bit); break; } }
                }
            }
            if (!anybit(b, nwd[Y])) goto out;
        } else if (cnt == 2) {
            starsum_dom(v, s);
            for (int side = 0; side < 2; side++) {
                int A = op[side], Bc = op[1 - side];
                u64 *bB = dlev[depth + 1][Bc], *bA = dlev[depth + 1][A];
                int nB = popc(bB, nwd[Bc]);
                int hs = 1;
                while (hs < 2 * nB + 2) hs <<= 1;
                if (hs > hcap) { free(htab); htab = malloc(sizeof(u64) * (size_t)hs); hcap = hs; }
                memset(htab, 0, sizeof(u64) * (size_t)hs);
                u64 *kb = (CC[Bc] == v) ? kp[Bc] : km[Bc];
                for (int w = 0; w < nwd[Bc]; w++) {
                    u64 x = bB[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        u64 key = kb[64 * w + bit] + 1ULL;
                        int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                        while (htab[pos] && htab[pos] != key) pos = (pos + 1) & (hs - 1);
                        htab[pos] = key;
                    }
                }
                for (int w = 0; w < nwd[A]; w++) {
                    u64 x = bA[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        const signed char *u = dom[A] + 34 * (64 * w + bit);
                        int needv[17], okr = 1;
                        for (int i = 0; i < 17; i++) {
                            int c = (CC[A] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                            needv[i] = -s[i] - c;
                            if (needv[i] < -4 || needv[i] > 4) { okr = 0; break; }
                        }
                        int found = 0;
                        if (okr) {
                            u64 key = ckey(NULL, 0, needv) + 1ULL;
                            int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                            while (htab[pos]) { if (htab[pos] == key) { found = 1; break; } pos = (pos + 1) & (hs - 1); }
                        }
                        if (!found) bA[w] &= ~(1ULL << bit);
                    }
                }
                if (!anybit(bA, nwd[A])) goto out;
            }
        }
    }
    res = rec_dom(depth + 1, R2);
out:
    assigned[X] = 0;
    return res;
}
static int rec_dom(int depth, double R[15][15]) {
    nodes++;
    if (nodes > cap) return 2;
    int best = -1, bc = 1 << 30;
    for (int Y = 0; Y < NC; Y++) if (!assigned[Y]) { int c = popc(dlev[depth][Y], nwd[Y]); if (c < bc) { bc = c; best = Y; } }
    if (best < 0) { for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0; return 1; }
    /* copy alive sets of assigned cells are irrelevant; unassigned ones are filtered into depth+1 by try_dom */
    u64 *b = dlev[depth][best];
    for (int w = 0; w < nwd[best]; w++) {
        u64 x = b[w];
        while (x) { int bit = __builtin_ctzll(x); x &= x - 1; int r = try_dom(depth, best, 64 * w + bit, R); if (r) return r; }
    }
    return 0;
}
static int try_bits(int k, int t, double R[15][15]) {
    if (!csign_ok(dom[k] + 34 * t, cfix[k], &cfix[k + 1])) return 0;
    const double *pv = proj[k] + 30 * t;
    double R2[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R2[i][j] = R[i][j] - pv[i] * pv[j] - pv[15 + i] * pv[15 + j];
    for (int i = 0; i < 15; i++) if (R2[i][i] < -1e-7) return 0;
    if (!psd_ok(R2)) return 0;
    chosen[k] = t;
    int s[17];
    for (int v = 0; v < 7; v++) if (lastk[v] == k) { starsum(v, k, s); for (int i = 0; i < 17; i++) if (s[i]) return 0; }
    if (!rowdone[k][t]) compat_row(k, t);
    for (int Y = k + 1; Y < NC; Y++) {
        const u64 *cm = compat[k][Y] + (size_t)t * nwd[Y];
        u64 *src = levb[k][Y], *dst = levb[k + 1][Y], acc = 0;
        for (int i = 0; i < nwd[Y]; i++) { dst[i] = src[i] & cm[i]; acc |= dst[i]; }
        if (!acc) return 0;
    }
    for (int v = 0; v < 7; v++) {
        int op[3], cnt = 0;
        for (int Y = k + 1; Y < NC; Y++) if (CC[Y] == v || CD[Y] == v) { if (cnt < 3) op[cnt] = Y; cnt++; }
        if (cnt == 1) {
            int Y = op[0];
            starsum(v, k, s);
            u64 *b = levb[k + 1][Y];
            for (int w = 0; w < nwd[Y]; w++) {
                u64 x = b[w];
                while (x) {
                    int bit = __builtin_ctzll(x); x &= x - 1;
                    const signed char *u = dom[Y] + 34 * (64 * w + bit);
                    for (int i = 0; i < 17; i++) { int c = (CC[Y] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i]; if (c != -s[i]) { b[w] &= ~(1ULL << bit); break; } }
                }
            }
            if (!anybit(b, nwd[Y])) return 0;
        } else if (cnt == 2) {
            starsum(v, k, s);
            for (int side = 0; side < 2; side++) {
                int A = op[side], Bc = op[1 - side];
                u64 *bB = levb[k + 1][Bc], *bA = levb[k + 1][A];
                int nB = 0;
                for (int w = 0; w < nwd[Bc]; w++) nB += __builtin_popcountll(bB[w]);
                int hs = 1;
                while (hs < 2 * nB + 2) hs <<= 1;
                if (hs > hcap) { free(htab); htab = malloc(sizeof(u64) * (size_t)hs); hcap = hs; }
                memset(htab, 0, sizeof(u64) * (size_t)hs);
                u64 *kb = (CC[Bc] == v) ? kp[Bc] : km[Bc];
                for (int w = 0; w < nwd[Bc]; w++) {
                    u64 x = bB[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        u64 key = kb[64 * w + bit] + 1ULL;
                        int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                        while (htab[pos] && htab[pos] != key) pos = (pos + 1) & (hs - 1);
                        htab[pos] = key;
                    }
                }
                for (int w = 0; w < nwd[A]; w++) {
                    u64 x = bA[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        const signed char *u = dom[A] + 34 * (64 * w + bit);
                        int needv[17], okr = 1;
                        for (int i = 0; i < 17; i++) {
                            int c = (CC[A] == v) ? u[i] + u[17 + i] : u[i] - u[17 + i];
                            needv[i] = -s[i] - c;
                            if (needv[i] < -4 || needv[i] > 4) { okr = 0; break; }
                        }
                        int found = 0;
                        if (okr) {
                            u64 key = ckey(NULL, 0, needv) + 1ULL;
                            int pos = (int)(((key * 0x9E3779B97F4A7C15ULL) >> 20) & (u64)(hs - 1));
                            while (htab[pos]) { if (htab[pos] == key) { found = 1; break; } pos = (pos + 1) & (hs - 1); }
                        }
                        if (!found) bA[w] &= ~(1ULL << bit);
                    }
                }
                if (!anybit(bA, nwd[A])) return 0;
            }
        } else if (cnt == 3) {
            starsum(v, k, s);
            for (int pr = 0; pr < 3; pr++) {
                if (!has3[v][pr]) continue;
                int A = L3v[v][pr];
                u64 *bA = levb[k + 1][A];
                for (int w = 0; w < nwd[A]; w++) {
                    u64 x = bA[w];
                    while (x) {
                        int bit = __builtin_ctzll(x); x &= x - 1;
                        int ca[17], need[17], okr = 1;
                        contrib(A, 64 * w + bit, v, ca);
                        for (int i = 0; i < 17; i++) { need[i] = -s[i] - ca[i]; if (need[i] < -16 || need[i] > 15) okr = 0; }
                        int found = 0;
                        if (okr) { u64 lo, hi; key85(need, &lo, &hi); found = h3_has(H3[v][pr], H3n[v][pr], lo, hi); }
                        if (!found) { bA[w] &= ~(1ULL << bit); st_sum3++; }
                    }
                }
                if (!anybit(bA, nwd[A])) return 0;
            }
        }
    }
    return rec_bits(k + 1, R2);
}
static int rec_bits(int k, double R[15][15]) {
    nodes++;
    if (nodes > cap) return 2;
    if (k == NC) { for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) if (fabs(R[i][j]) > 1e-6) return 0; return 1; }
    u64 *b = levb[k][k];
    for (int w = 0; w < nwd[k]; w++) {
        u64 x = b[w];
        while (x) { int bit = __builtin_ctzll(x); x &= x - 1; int r = try_bits(k, 64 * w + bit, R); if (r) return r; }
    }
    return 0;
}

static void prep_common(void) {
    for (int X = 0; X < NC; X++) {
        int n = ndom[X];
        proj[X] = malloc(sizeof(double) * 30 * (size_t)n);
        for (int t = 0; t < n; t++) for (int s2 = 0; s2 < 2; s2++) for (int j = 0; j < 15; j++) {
            double a = 0;
            for (int i = 0; i < 17; i++) a += dom[X][34 * t + 17 * s2 + i] * Bw[i * 15 + j];
            proj[X][30 * t + 15 * s2 + j] = a;
        }
    }
}
static void free_common_dom(void) { for (int X = 0; X < NC; X++) { free(dom[X]); dom[X] = NULL; free(ipair[X]); ipair[X] = NULL; } }
static void free_common(void) { for (int X = 0; X < NC; X++) { free(proj[X]); free(dom[X]); dom[X] = NULL; free(ipair[X]); ipair[X] = NULL; } }

static int solve_plain(long long capv, long long *nodes_out) {
    for (int X = 0; X < NC; X++) {
        for (int L = 0; L <= NC; L++) levi[L][X] = malloc(sizeof(int) * (size_t)ndom[X]);
        for (int t = 0; t < ndom[X]; t++) levi[0][X][t] = t;
        levn[0][X] = ndom[X];
    }
    nodes = 0; cap = capv;
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    int r = rec_plain(0, R);
    for (int X = 0; X < NC; X++) for (int L = 0; L <= NC; L++) free(levi[L][X]);
    *nodes_out = nodes;
    return r;
}
static int solve_bits(long long capv, long long *nodes_out) {
    for (int X = 0; X < NC; X++) {
        nwd[X] = (ndom[X] + 63) / 64;
        kp[X] = malloc(sizeof(u64) * (size_t)ndom[X]); km[X] = malloc(sizeof(u64) * (size_t)ndom[X]);
        for (int t = 0; t < ndom[X]; t++) { kp[X][t] = ckey(dom[X] + 34 * t, 1, NULL); km[X][t] = ckey(dom[X] + 34 * t, 0, NULL); }
        for (int L = 0; L <= NC; L++) levb[L][X] = calloc((size_t)nwd[X], sizeof(u64));
        for (int t = 0; t < ndom[X]; t++) levb[0][X][t >> 6] |= 1ULL << (t & 63);
    }
    for (int X = 0; X < NC; X++) {
        rowdone[X] = calloc((size_t)(ndom[X] ? ndom[X] : 1), 1);
        for (int Y = X + 1; Y < NC; Y++) compat[X][Y] = calloc((size_t)ndom[X] * nwd[Y], sizeof(u64));
    }
    for (int v = 0; v < 7; v++) {
        int cnt = 0, cells[6];
        for (int X = 0; X < NC; X++) if (CC[X] == v || CD[X] == v) cells[cnt++] = X;
        L3v[v][0] = cells[3]; L3v[v][1] = cells[4]; L3v[v][2] = cells[5];
        for (int pr = 0; pr < 3; pr++) {                   /* pr = index of the EXCLUDED cell; the other two form the pair */
            int P = L3v[v][pr == 0 ? 1 : 0], Q = L3v[v][pr == 2 ? 1 : 2];
            has3[v][pr] = 0; H3[v][pr] = NULL;
            if ((long long)ndom[P] * ndom[Q] > look3_cap) continue;
            long long n = 3 * (long long)ndom[P] * ndom[Q] + 16;
            u64 *T = calloc((size_t)n * 2, sizeof(u64));
            for (int a = 0; a < ndom[P]; a++) {
                const signed char *wa = dom[P] + 34 * a;
                int ca[17]; contrib(P, a, v, ca);
                for (int b = 0; b < ndom[Q]; b++) {
                    const signed char *wb = dom[Q] + 34 * b;
                    int ok = 1;
                    for (int s2 = 0; s2 < 2 && ok; s2++) for (int tt = 0; tt < 2 && ok; tt++) {
                        int g = rrt[P][Q][s2][tt];
                        for (int i = 0; i < 17; i++) g += wa[17 * s2 + i] * wb[17 * tt + i];
                        if (!(g == 0 || g == 4 || g == -4)) ok = 0;
                    }
                    if (!ok) continue;
                    int cb[17], sm[17]; contrib(Q, b, v, cb);
                    for (int i = 0; i < 17; i++) sm[i] = ca[i] + cb[i];
                    u64 lo, hi; key85(sm, &lo, &hi); h3_insert(T, n, lo, hi);
                }
            }
            H3[v][pr] = T; H3n[v][pr] = n; has3[v][pr] = 1;
        }
    }
    nodes = 0; cap = capv;
    double R[15][15];
    for (int i = 0; i < 15; i++) for (int j = 0; j < 15; j++) R[i][j] = (i == j) ? 28.0 : 0.0;
    int r;
    if (use_ac && !ac_full()) { r = 0; nodes = 1; st_acwipe++; }
    else if (use_dom) {
        for (int X = 0; X < NC; X++) for (int t = 0; t < ndom[X]; t++) if (!rowdone[X][t]) compat_row(X, t);
        for (int X = 0; X < NC; X++) for (int Y = 0; Y < X; Y++) {
            cmT[X][Y] = calloc((size_t)ndom[X] * nwd[Y], sizeof(u64));
            for (int q = 0; q < ndom[Y]; q++) {
                const u64 *row = compat[Y][X] + (size_t)q * nwd[X];
                for (int t = 0; t < ndom[X]; t++) if (row[t >> 6] >> (t & 63) & 1) cmT[X][Y][(size_t)t * nwd[Y] + (q >> 6)] |= 1ULL << (q & 63);
            }
        }
        for (int L = 0; L <= NC; L++) for (int X = 0; X < NC; X++) {
            dlev[L][X] = calloc((size_t)nwd[X], sizeof(u64));
        }
        if (use_zsac && use_z4 && use_zfilter && nZM > 0 && !zmOverflow) {
            /* per-code restriction: entries with the forced signature, then AC, then search (complete: each solution lies in one code) */
            u64 *Eb[NC][64];
            for (int X = 0; X < NC; X++) for (int sg = 0; sg < 64; sg++) Eb[X][sg] = NULL;
            for (int X = 0; X < NC; X++) for (int t2 = 0; t2 < ndom[X]; t2++) {
                int sg = sigE[X][t2];
                if (!Eb[X][sg]) Eb[X][sg] = calloc((size_t)nwd[X], sizeof(u64));
                Eb[X][sg][t2 >> 6] |= 1ULL << (t2 & 63);
            }
            u64 *al[NC];
            for (int X = 0; X < NC; X++) al[X] = malloc(sizeof(u64) * (size_t)nwd[X]);
            r = 0;
            long long tot = 0;
            for (int m = 0; m < nZM && r != 1; m++) {
                int ok = 1;
                for (int X = 0; X < NC && ok; X++) {
                    const u64 *eb = Eb[X][zT[(size_t)m * NC + X]];
                    if (!eb) { ok = 0; break; }
                    u64 acc = 0;
                    for (int i = 0; i < nwd[X]; i++) { al[X][i] = levb[0][X][i] & eb[i]; acc |= al[X][i]; }
                    if (!acc) ok = 0;
                }
                st_zsac_codes++;
                if (!ok || !ac_bits(al)) { st_zsac_acdead++; continue; }
                st_zsac_search++;
                for (int X = 0; X < NC; X++) { memcpy(dlev[0][X], al[X], sizeof(u64) * (size_t)nwd[X]); assigned[X] = 0; }
                nodes = 0;
                int rr2 = rec_dom(0, R);
                tot += nodes;
                if (rr2 == 1) r = 1;
                else if (rr2 == 2 && r == 0) r = 2;
            }
            nodes = tot + 1;
            for (int X = 0; X < NC; X++) { free(al[X]); for (int sg = 0; sg < 64; sg++) free(Eb[X][sg]); }
        } else {
            for (int X = 0; X < NC; X++) { memcpy(dlev[0][X], levb[0][X], sizeof(u64) * (size_t)nwd[X]); assigned[X] = 0; }
            r = rec_dom(0, R);
        }
        for (int L = 0; L <= NC; L++) for (int X = 0; X < NC; X++) free(dlev[L][X]);
        for (int X = 0; X < NC; X++) for (int Y = 0; Y < X; Y++) { free(cmT[X][Y]); cmT[X][Y] = NULL; }
    }
    else r = rec_bits(0, R);
    for (int v = 0; v < 7; v++) for (int pr = 0; pr < 3; pr++) if (has3[v][pr]) { free(H3[v][pr]); has3[v][pr] = 0; }
    for (int X = 0; X < NC; X++) {
        free(kp[X]); free(km[X]);
        for (int L = 0; L <= NC; L++) free(levb[L][X]);
        for (int Y = X + 1; Y < NC; Y++) free(compat[X][Y]);
        free(rowdone[X]);
    }
    *nodes_out = nodes;
    return r;
}

/* ---------------- S_7 x class-map canonical form (dedupe) ---------------- */
static int dedupe = 0, nmaps = 1, cmap[8][64], ncls, clsn[64], clsidx[64][17];
static signed char *PT;            /* 5040 x 128 */
static int nperm = 0;
static u64 *hset; static long long hsize = 0, hcount = 0, st_dup = 0;
static void gen_perms(void) {
    int p[7] = {0, 1, 2, 3, 4, 5, 6};
    PT = malloc(5040 * 128);
    int c[7] = {0}, i = 0;
    for (;;) {
        for (int x = 0; x < 128; x++) {
            int y = 0;
            for (int b = 0; b < 7; b++) if (x >> b & 1) y |= 1 << p[b];
            if (y & 1) y ^= 127;
            PT[nperm * 128 + x] = (signed char)y;
        }
        nperm++;
        while (i < 7 && c[i] >= i) { c[i] = 0; i++; }
        if (i >= 7) break;
        if (i % 2 == 0) { int t = p[0]; p[0] = p[i]; p[i] = t; } else { int t = p[c[i]]; p[c[i]] = p[i]; p[i] = t; }
        c[i]++; i = 0;
    }
}
static void canon_form(const int *cols, unsigned char *best) {
    unsigned char cur[17];
    int first = 1;
    for (int pi = 0; pi < nperm; pi++) {
        const signed char *T = PT + pi * 128;
        for (int mi = 0; mi < nmaps; mi++) {
            int pos = 0;
            for (int a = 0; a < ncls; a++) {
                int src = cmap[mi][a];
                unsigned char v[17];
                int n = clsn[src];
                for (int j = 0; j < n; j++) v[j] = (unsigned char)T[cols[clsidx[src][j]]];
                for (int x = 1; x < n; x++) { unsigned char t = v[x]; int y = x - 1; while (y >= 0 && v[y] > t) { v[y + 1] = v[y]; y--; } v[y + 1] = t; }
                for (int j = 0; j < n; j++) cur[pos++] = v[j];
            }
            if (first || memcmp(cur, best, 17) < 0) { memcpy(best, cur, 17); first = 0; }
        }
    }
}
static int hs_insert(const unsigned char *k17) {
    u64 k[3] = {0, 0, 0};
    memcpy(k, k17, 17);
    k[2] |= 1ULL << 63;
    if (hcount * 2 >= hsize) {
        long long ns = hsize ? hsize * 2 : (1 << 20);
        u64 *nh = calloc((size_t)ns * 3, sizeof(u64));
        for (long long i = 0; i < hsize; i++) if (hset[3 * i + 2]) {
            u64 h = (hset[3 * i] * 0x9E3779B97F4A7C15ULL) ^ (hset[3 * i + 1] * 0xC2B2AE3D27D4EB4FULL) ^ hset[3 * i + 2];
            long long q = (long long)(h % (u64)ns);
            while (nh[3 * q + 2]) q = (q + 1) % ns;
            memcpy(nh + 3 * q, hset + 3 * i, 24);
        }
        free(hset); hset = nh; hsize = ns;
    }
    u64 h = (k[0] * 0x9E3779B97F4A7C15ULL) ^ (k[1] * 0xC2B2AE3D27D4EB4FULL) ^ k[2];
    long long q = (long long)(h % (u64)hsize);
    while (hset[3 * q + 2]) {
        if (hset[3 * q] == k[0] && hset[3 * q + 1] == k[1] && hset[3 * q + 2] == k[2]) return 0;
        q = (q + 1) % hsize;
    }
    memcpy(hset + 3 * q, k, 24); hcount++;
    return 1;
}

/* ---------------- stage A ---------------- */
typedef struct { int n; int idx[17]; } Part;
static unsigned rows[7];
static int nrows;
/* Knuth estimator of the stage-A tree (mode "knuth NW seed [solve]"): random root-to-leaf walks, one child drawn uniformly
   per level by reservoir sampling; the product of the branching factors is an unbiased estimate of the number of leaves. */
static int knuth_mode = 0, knuth_solve = 0, knuth_pz4 = 0;
static double kn_pzkill[8], kn_pzsec[8];
static int kn_firstfail;
static long long kn_cnt[8]; static unsigned kn_m[8]; static Part kn_parts[8][64]; static int kn_np[8];

static int m1_ok(void) {
    for (int c = 0; c < nrows; c++)
        for (int d = c + 1; d < nrows; d++) {
            unsigned S = rows[c] ^ rows[d];
            int *L = amask[S], n = namask[S], found = 0;
            for (int t = 0; t < n && !found; t++) {
                int good = 1;
                for (int e = 1; e < nrows && good; e++) {
                    int par = pc((unsigned)L[t] & rows[e]) & 1, tg = (c == 0) ? 1 : (e == c);
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
        for (int W = w; W <= w + rem && !ok; W++) { int k = W < 7 - W ? W : 7 - W; if (k >= lo && k <= peK[e]) ok = 1; }
        if (!ok) return 0;
    }
    return 1;
}

static void leaf(void) {
    if (maxtime && time(NULL) - t_start >= maxtime) { stop_now = 1; return; }
    int cols[17];
    for (int i = 0; i < 17; i++) {
        int c = 0;
        for (int r = 0; r < 7; r++) if (rows[r] >> i & 1) c |= 1 << r;
        if (!c) return;
        cols[i] = c;
    }
    st_leaves++;
    if (samp < 1.0 && urand() > samp) return;
    if (dedupe) {
        unsigned char cf[17];
        canon_form(cols, cf);
        if (!hs_insert(cf)) { st_dup++; return; }
    }
    if (countonly) {
        st_proc++;
        if (time(NULL) - t_last >= 60) {
            t_last = time(NULL);
            printf("  count: top %lld/%lld dup %lld distinct %lld leaves %lld (%lds)%c", top_counter, top_total, st_dup, st_proc, st_leaves,
                   (long)(time(NULL) - t_start), 10);
            fflush(stdout);
        }
        return;
    }
    st_proc++;
    clock_t c0 = clock();
    int bok = build(rows, cols);
    tm_build += (double)(clock() - c0) / CLOCKS_PER_SEC;
    if (bok && use_z4 && !z4_ok(rows)) {
        st_z4++;
        if (0) {}
        if (fdbg && st_proc <= debugN) { for (int i = 0; i < 17; i++) fprintf(fdbg, "%02x", cols[i]); fprintf(fdbg, " Z4DEAD 0%c", 10); }
        free_common_dom();
        goto progress;
    }
    if (bok && use_z4 && use_zfilter) zfilter();
    if (!bok) {
        st_empty++;
        if (fdbg && st_proc <= debugN) { for (int i = 0; i < 17; i++) fprintf(fdbg, "%02x", cols[i]); fprintf(fdbg, " EMPTY 0\n"); }
        goto progress;
    }
    c0 = clock();
    prep_common();
    tm_prep += (double)(clock() - c0) / CLOCKS_PER_SEC;
    c0 = clock();
    long long nd;
    int r = solve_plain(2000, &nd);
    if (r == 2) {
        r = solve_bits(heavycap, &nd);
        if (r == 0 && nd > 100000) {
            st_heavy++;
            fprintf(fsurv, "HEAVY-EXCLUDED %lld ", nd);
            for (int i = 0; i < 17; i++) fprintf(fsurv, "%02x", cols[i]);
            fprintf(fsurv, "\n"); fflush(fsurv);
        }
    }
    tm_solve += (double)(clock() - c0) / CLOCKS_PER_SEC;
    st_nodes += nd;
    if (fdbg && st_proc <= debugN) { for (int i = 0; i < 17; i++) fprintf(fdbg, "%02x", cols[i]); fprintf(fdbg, " %s %lld\n", r == 0 ? "EXCLUDED" : (r == 1 ? "FOUND" : "LIMIT"), nd); }
    if (r == 0) { if (nd == 1) st_root++; else st_deeper++; }
    else {
        if (r == 1) st_found++; else st_limit++;
        fprintf(fsurv, "%s %lld ", r == 1 ? "FOUND" : "LIMIT", nd);
        for (int i = 0; i < 17; i++) fprintf(fsurv, "%02x", cols[i]);
        fprintf(fsurv, "\n");
        if (r == 1)
            for (int X = 0; X < NC; X++) {
                fprintf(fsurv, "  %d %d", CC[X], CD[X]);
                for (int i = 0; i < 34; i++) fprintf(fsurv, " %d", dom[X][34 * chosen[X] + i]);
                fprintf(fsurv, "\n");
            }
        fflush(fsurv);
    }
    free_common();
progress:
    if (time(NULL) - t_last >= 120) {
        t_last = time(NULL);
        printf("  zsac codes %lld acdead %lld searched %lld | zf %lld->%lld (M %lld) z4 %lld dup %lld | leaves %lld processed %lld empty %lld root %lld deeper %lld heavy %lld found %lld limit %lld m1_prunes %lld (%lds) t_build %.0f t_prep %.0f t_solve %.0f\n",
               st_zsac_codes, st_zsac_acdead, st_zsac_search, st_zfilt_in, st_zfilt_out, st_zM, st_z4, st_dup, st_leaves, st_proc, st_empty, st_root, st_deeper, st_heavy, st_found, st_limit, m1_prunes, (long)(time(NULL) - t_start),
               tm_build, tm_prep, tm_solve);
        fflush(stdout);
    }
}

static int GOOD(int w) { return w == 2 || w == 6 || w == 10; }
static void dfs(Part *parts, int np);
static int selv[64];
static void choose(Part *parts, int np, int a, int wt, unsigned forced) {
    if (wt > 10 || stop_now) return;
    if (a == np) {
        if (!GOOD(wt)) return;
        int c = nrows, w2 = c >= 2 ? pc(rows[1]) : 2;
        int wmin = (c <= 2) ? w2 : (w2 > pc(rows[c - 1]) ? w2 : pc(rows[c - 1]));
        if (c == 1) wmin = 2;
        if (wt < wmin) return;
        unsigned m = 0;
        for (int b = 0; b < np; b++) for (int j = 0; j < selv[b]; j++) m |= 1u << parts[b].idx[j];
        if ((pc(m & (unsigned)O1) & 1) || (pc(m & (unsigned)O2) & 1)) return;
        if (!usable(m)) return;
        for (int r = 1; r < c; r++) { int d = pc(m ^ rows[r]); if (!GOOD(d) || (c >= 2 && d < w2) || !usable(m ^ rows[r])) return; }
        if (c == 1) {                          /* top-level branch: work split across machines */
            if (toponly) { top_total++; return; }
            long long tc = top_counter++;
            if (tc % part_m != part_r) return;
        }
        if (probe_q < 1.0 && c >= 2 && urand() > probe_q) return;     /* random probing: keep each branch with prob q */
        rows[nrows++] = m;
        int ok = m1_ok();
        if (!ok) m1_prunes++;
        if (ok && kfeasible()) {
            Part np2[64];
            int n2 = 0;
            for (int b = 0; b < np; b++) {
                if (selv[b]) { np2[n2].n = selv[b]; memcpy(np2[n2].idx, parts[b].idx, sizeof(int) * selv[b]); n2++; }
                if (selv[b] < parts[b].n) { np2[n2].n = parts[b].n - selv[b]; memcpy(np2[n2].idx, parts[b].idx + selv[b], sizeof(int) * (parts[b].n - selv[b])); n2++; }
            }
            if (knuth_mode) {
                int L = nrows - 1;
                kn_cnt[L]++;
                if (urand() * (double)kn_cnt[L] < 1.0) { kn_m[L] = m; kn_np[L] = n2; memcpy(kn_parts[L], np2, sizeof(Part) * (size_t)n2); }
            } else dfs(np2, n2);
        }
        nrows--;
        return;
    }
    if (forced >> a & 1) { selv[a] = parts[a].n; choose(parts, np, a + 1, wt + parts[a].n, forced); selv[a] = 0; return; }
    for (int s = 0; s <= parts[a].n; s++) { selv[a] = s; choose(parts, np, a + 1, wt + s, forced); }
    selv[a] = 0;
}
static void dfs(Part *parts, int np) {
    if (nrows == 7) { leaf(); return; }
    unsigned forced = 0;
    if (nrows == 6) {
        unsigned cov = 0;
        for (int r = 0; r < 6; r++) cov |= rows[r];
        for (int a = 0; a < np; a++) if (!(cov >> parts[a].idx[0] & 1)) forced |= 1u << a;
    }
    int save[64];
    memcpy(save, selv, sizeof save); memset(selv, 0, sizeof selv);
    choose(parts, np, 0, 0, forced);
    memcpy(selv, save, sizeof save);
}

/* ---- partial Z4 test at stage-A depth k (sound: M1 and the signatures restricted to rows 1..k, no L1 refinement) ---- */
typedef struct { unsigned x, am; unsigned char fl; } PTrip;
static PTrip **ptr_c; static int *ptr_n;
static void ptrip_build(unsigned S) {
    if (ptr_c[S]) return;
    usable(S);
    int n = ncvec[S];
    const signed char *W = cvec[S];
    int cap = 256, m = 0;
    PTrip *T = malloc(sizeof(PTrip) * (size_t)cap);
    /* hash of distinct triples */
    int hc = 1 << 12;
    u64 *H = calloc((size_t)hc, sizeof(u64));
    for (int p = 0; p < n; p++) {
        const signed char *wp = W + 17 * p;
        unsigned x = 0;
        for (int i = 0; i < 17; i++) { int q = wp[i], xi; if (q & 1) xi = ((((q - 1) % 4) + 4) % 4) / 2; else xi = (((q % 4) + 4) % 4) / 2; if (xi) x |= 1u << i; }
        for (int q = 0; q < n; q++) {
            if (q == p) continue;
            const signed char *wq = W + 17 * q;
            int g = 0;
            for (int i = 0; i < 17; i++) g += wp[i] * wq[i];
            if (!(g == 0 || g == 4 || g == -4)) continue;
            int a[17], b[17]; unsigned am = 0;
            for (int i = 0; i < 17; i++) { a[i] = (wp[i] + wq[i]) / 2; b[i] = (wp[i] - wq[i]) / 2; if (a[i] & 1) am |= 1u << i; }
            unsigned char fl = 0;
            if (use_csign) { if (first_pos_nz(b) != 0) fl |= 1; if (first_pos_nz(a) != 0) fl |= 2; }
            else { if (first_pos(b)) fl |= 1; if (first_pos(a)) fl |= 2; }
            u64 key = ((u64)x << 20 | am) << 2 | fl, h = (key * 0x9E3779B97F4A7C15ULL) >> 20;
            int pos = (int)(h & (u64)(hc - 1)), dup = 0;
            while (H[pos]) { if (H[pos] == key + 1) { dup = 1; break; } pos = (pos + 1) & (hc - 1); }
            if (dup) continue;
            H[pos] = key + 1;
            if (m == cap) { cap *= 2; T = realloc(T, sizeof(PTrip) * (size_t)cap); }
            T[m].x = x; T[m].am = am; T[m].fl = fl; m++;
            if (2 * m > hc) {                           /* grow and rehash */
                int hc2 = hc * 2;
                u64 *H2 = calloc((size_t)hc2, sizeof(u64));
                for (int j = 0; j < m; j++) {
                    u64 k2 = ((u64)T[j].x << 20 | T[j].am) << 2 | T[j].fl, h2 = (k2 * 0x9E3779B97F4A7C15ULL) >> 20;
                    int ps = (int)(h2 & (u64)(hc2 - 1));
                    while (H2[ps]) ps = (ps + 1) & (hc2 - 1);
                    H2[ps] = k2 + 1;
                }
                free(H); H = H2; hc = hc2;
            }
        }
    }
    free(H);
    ptr_c[S] = T; ptr_n[S] = m;
}
static u64 pSig[7][7];
static int pM[7], pk;
static int pJ[7][7], pkap[7][7];
static int pz4_rec(int c) {
    if (c > pk) return 1;
    u64 S = pSig[0][c];
    while (S) {
        int m = __builtin_ctzll(S); S &= S - 1;
        int ok = 1;
        for (int e = 1; e < c && ok; e++) if ((((m >> (e - 1)) & 1) ^ ((pM[e] >> (c - 1)) & 1)) != pJ[c][e]) ok = 0;
        for (int e = 1; e < c && ok; e++) { int v = pM[e] ^ m ^ pkap[e][c]; if (!(pSig[e][c] >> v & 1)) ok = 0; }
        if (!ok) continue;
        pM[c] = m;
        if (pz4_rec(c + 1)) return 1;
    }
    return 0;
}
static long long st_pz4_calls = 0, st_pz4_kill = 0;
static int pz4_ok(int k) {                 /* rows[0..k] set; tests the cells among vertices 0..k */
    st_pz4_calls++;
    unsigned bzk[7];
    for (int c = 1; c <= k; c++) bzk[c] = 1u | (1u << c) | (rows[c] << 7);
    for (int c = 1; c <= k; c++) for (int e = 1; e <= k; e++) pJ[c][e] = (pc(bzk[c] & bzk[e]) / 2) & 1;
    for (int c = 1; c <= k; c++) for (int d = 1; d <= k; d++) {
        int v = 0;
        for (int e = 1; e <= k; e++) if (pc(bzk[c] & bzk[d] & bzk[e]) & 1) v |= 1 << (e - 1);
        pkap[c][d] = v;
    }
    for (int c = 0; c <= k; c++) for (int d = c + 1; d <= k; d++) {
        unsigned S = rows[c] ^ rows[d];
        ptrip_build(S);
        const PTrip *T = ptr_c[S]; int n = ptr_n[S];
        unsigned char need = (c == 0) ? (d == 1 ? 3 : 1) : 0;
        int tgv = 0;                                  /* required parities of a.rows[e], e = 1..k */
        for (int e = 1; e <= k; e++) { int tg = (c == 0) ? 1 : (e == c); if (tg) tgv |= 1 << (e - 1); }
        u64 sig = 0;
        for (int t = 0; t < n; t++) {
            if ((T[t].fl & need) != need) continue;
            int par = 0, sg = 0;
            for (int e = 1; e <= k; e++) {
                if (pc(T[t].am & rows[e]) & 1) par |= 1 << (e - 1);
                if (pc(T[t].x & rows[e]) & 1) sg |= 1 << (e - 1);
            }
            if (par != tgv) continue;
            sig |= 1ULL << sg;
        }
        if (!sig) { st_pz4_kill++; return 0; }
        pSig[c][d] = sig;
    }
    pk = k;
    int r = pz4_rec(1);
    if (!r) st_pz4_kill++;
    return r;
}

/* one Knuth walk from the root; returns the weight of the reached leaf (0 on a dead end); *cls: 0 dead end or uncovered,
   1 empty/z4-dead, 2 excluded by the solver, 3 found/limit (only with knuth_solve); *secs: solver seconds */
static double knuth_walk(Part *parts0, int np0, int *cls, double *secs, double *wdepth) {
    Part cur[64]; int ncur = np0;
    memcpy(cur, parts0, sizeof(Part) * (size_t)np0);
    double w = 1.0;
    rows[0] = 0; nrows = 1;
    *cls = 0; *secs = 0; kn_firstfail = 0;
    while (nrows < 7) {
        int L = nrows;
        kn_cnt[L] = 0;
        unsigned forced = 0;
        if (nrows == 6) {
            unsigned cov = 0;
            for (int r = 0; r < 6; r++) cov |= rows[r];
            for (int a = 0; a < ncur; a++) if (!(cov >> cur[a].idx[0] & 1)) forced |= 1u << a;
        }
        memset(selv, 0, sizeof selv);
        choose(cur, ncur, 0, 0, forced);
        if (kn_cnt[L] == 0) return 0.0;
        w *= (double)kn_cnt[L];
        wdepth[L] += w;
        rows[nrows++] = kn_m[L];
        ncur = kn_np[L];
        memcpy(cur, kn_parts[L], sizeof(Part) * (size_t)ncur);
        if (knuth_pz4 && L >= 2) {
            clock_t c1 = clock();
            int ok = pz4_ok(L);
            kn_pzsec[L] += w * (double)(clock() - c1) / CLOCKS_PER_SEC;
            if (!ok) { kn_pzkill[L] += w; if (!kn_firstfail) kn_firstfail = L; }
        }
    }
    int cols[17];
    for (int i = 0; i < 17; i++) {
        int c = 0;
        for (int r = 0; r < 7; r++) if (rows[r] >> i & 1) c |= 1 << r;
        if (!c) return 0.0;
        cols[i] = c;
    }
    *cls = 1;
    if (!knuth_solve) return w;
    clock_t c0 = clock();
    if (!build(rows, cols)) { *secs = (double)(clock() - c0) / CLOCKS_PER_SEC; return w; }
    if (use_z4 && !z4_ok(rows)) { free_common_dom(); *secs = (double)(clock() - c0) / CLOCKS_PER_SEC; return w; }
    if (use_z4 && use_zfilter) zfilter();
    prep_common();
    long long nd; int r = solve_plain(2000, &nd);
    if (r == 2) r = solve_bits(heavycap, &nd);
    free_common();
    *secs = (double)(clock() - c0) / CLOCKS_PER_SEC;
    *cls = (r == 0) ? 2 : 3;
    return w;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    char fn[1024];
    for (int a = 3; a < argc; a++) {
        if (!strcmp(argv[a], "sample") && a + 2 < argc) { samp = atof(argv[a + 1]); rng_state ^= (u64)atoll(argv[a + 2]) * 0x9E3779B97F4A7C15ULL; if (!rng_state) rng_state = 1; a += 2; }
        else if (!strcmp(argv[a], "part") && a + 2 < argc) { part_r = atoi(argv[a + 1]); part_m = atoi(argv[a + 2]); a += 2; }
        else if (!strcmp(argv[a], "heavycap") && a + 1 < argc) { heavycap = atoll(argv[a + 1]); a += 1; }
        else if (!strcmp(argv[a], "dedupe")) { dedupe = 1; }
        else if (!strcmp(argv[a], "countonly")) { countonly = 1; }
        else if (!strcmp(argv[a], "noac")) { use_ac = 0; }
        else if (!strcmp(argv[a], "nocsign")) { use_csign = 0; }
        else if (!strcmp(argv[a], "maxtime") && a + 1 < argc) { maxtime = atol(argv[a + 1]); a += 1; }
        else if (!strcmp(argv[a], "nodom")) { use_dom = 0; }
        else if (!strcmp(argv[a], "noz4")) { use_z4 = 0; }
        else if (!strcmp(argv[a], "nozfilter")) { use_zfilter = 0; }
        else if (!strcmp(argv[a], "nozsac")) { use_zsac = 0; }
        else if (!strcmp(argv[a], "probe") && a + 2 < argc) { probe_q = atof(argv[a + 1]); rng_state ^= (u64)atoll(argv[a + 2]) * 0xD1B54A32D192ED03ULL; if (!rng_state) rng_state = 7; a += 2; }
        else if (!strcmp(argv[a], "knuth") && a + 2 < argc) { knuth_mode = atoi(argv[a + 1]); rng_state ^= (u64)atoll(argv[a + 2]) * 0xA0761D6478BD642FULL; if (!rng_state) rng_state = 3; a += 2; }
        else if (!strcmp(argv[a], "solve")) { knuth_solve = 1; }
        else if (!strcmp(argv[a], "pz4")) { knuth_pz4 = 1; }
        else if (!strcmp(argv[a], "debug") && a + 2 < argc) { debugN = atoll(argv[a + 1]); fdbg = fopen(argv[a + 2], "w"); a += 2; }
    }
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
    for (int j = 0; j < 17 * 15; j++) if (fscanf(f, "%lf", &Bw[j]) != 1) { fprintf(stderr, "need Bw in input\n"); return 1; }
    ncls = ncl;
    for (int a = 0; a < ncl; a++) { clsn[a] = parts[a].n; memcpy(clsidx[a], parts[a].idx, sizeof(int) * parts[a].n); }
    nmaps = 1;
    for (int a = 0; a < ncl; a++) cmap[0][a] = a;
    {
        int nm;
        if (fscanf(f, "%d", &nm) == 1 && nm >= 1) {
            nmaps = nm;
            for (int mi = 0; mi < nm; mi++) for (int a = 0; a < ncl; a++) if (fscanf(f, "%d", &cmap[mi][a]) != 1) return 1;
        }
    }
    fclose(f);
    usab = calloc(1u << 17, 1); amask = calloc(1u << 17, sizeof(int *)); namask = calloc(1u << 17, sizeof(int));
    ptr_c = calloc(1u << 17, sizeof(PTrip *)); ptr_n = calloc(1u << 17, sizeof(int));
    cvec = calloc(1u << 17, sizeof(signed char *)); ncvec = calloc(1u << 17, sizeof(int));
    int k = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { CC[k] = c; CD[k] = d; k++; }
    for (int X = 0; X < NC; X++) for (int Y = 0; Y < NC; Y++) for (int s = 0; s < 2; s++) for (int t = 0; t < 2; t++) {
        int r1[7] = {0}, r2[7] = {0}, v = 0;
        r1[CC[X]] = 1; r1[CD[X]] = s ? -1 : 1; r2[CC[Y]] = 1; r2[CD[Y]] = t ? -1 : 1;
        for (int i = 0; i < 7; i++) v += r1[i] * r2[i];
        rrt[X][Y][s][t] = v;
    }
    for (int v = 0; v < 7; v++) lastk[v] = -1;
    for (int x = 0; x < NC; x++) { lastk[CC[x]] = x; lastk[CD[x]] = x; }
    for (int i = 0; i < 17; i++) if (Y1[i] == 0 && Y2[i] == 0) zc_mask |= 1u << i;
    cfix[0] = 0;
    snprintf(fn, sizeof fn, "%s.surv", argv[2]); fsurv = fopen(fn, "w");
    if (dedupe) gen_perms();
    for (int a = 3; a < argc; a++) if (!strcmp(argv[a], "leaves") && a + 1 < argc) {
        FILE *ft = fopen(argv[a + 1], "r");
        char line[128];
        while (fgets(line, sizeof line, ft)) {
            int cols[17];
            for (int i = 0; i < 17; i++) { unsigned v; sscanf(line + 2 * i, "%2x", &v); cols[i] = (int)v; }
            unsigned rw[7];
            for (int e = 0; e < 7; e++) { rw[e] = 0; for (int i = 0; i < 17; i++) if (cols[i] >> e & 1) rw[e] |= 1u << i; }
            for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) usable(rw[c] ^ rw[d]);
            clock_t c0 = clock();
            if (!build(rw, cols)) { printf("EMPTY 0 0.0s%c", 10); continue; }
            if (use_z4 && !z4_ok(rw)) { printf("Z4DEAD 0 0.0s%c", 10); free_common_dom(); continue; }
            if (use_z4 && use_zfilter) zfilter();
            prep_common();
            long long nd; int r = solve_plain(2000, &nd), tier = 1;
            if (r == 2) { r = solve_bits(heavycap, &nd); tier = 2; }
            printf("%s %lld tier%d %.1fs sum3prunes %lld%c", r == 0 ? "EXCLUDED" : (r == 1 ? "FOUND" : "LIMIT"), nd, tier,
                   (double)(clock() - c0) / CLOCKS_PER_SEC, st_sum3, 10);
            fflush(stdout);
            free_common();
        }
        fclose(ft);
        return 0;
    }
    for (int a = 3; a < argc; a++) if (!strcmp(argv[a], "canontest") && a + 1 < argc) {
        if (!dedupe) gen_perms();
        FILE *ft = fopen(argv[a + 1], "r");
        char line[128];
        while (fgets(line, sizeof line, ft)) {
            int cols[17];
            for (int i = 0; i < 17; i++) { unsigned v; sscanf(line + 2 * i, "%2x", &v); cols[i] = (int)v; }
            unsigned char cf[17];
            canon_form(cols, cf);
            for (int i = 0; i < 17; i++) printf("%d%s", cf[i], i < 16 ? " " : "");
            printf("\n");
        }
        fclose(ft);
        return 0;
    }
    if (knuth_mode) {
        double sw = 0, sw2 = 0, wz = 0, wz2 = 0, ws = 0, ws2 = 0, wsec = 0, wsec2 = 0, wdepth[8] = {0}, wpass = 0, wpass2 = 0;
        long long nleaf = 0, nsurv = 0;
        time_t t0 = time(NULL);
        for (int k = 1; k <= knuth_mode; k++) {
            int cls; double secs;
            double w = knuth_walk(parts, ncl, &cls, &secs, wdepth);
            sw += w; sw2 += w * w;
            if (cls >= 1) nleaf++;
            if (cls >= 1 && !kn_firstfail) { wpass += w; wpass2 += w * w; }
            double z = (cls >= 2) ? w : 0.0;       /* leaves reaching the integer solver (z4 survivors) */
            wz += z; wz2 += z * z;
            double sv = (cls == 3) ? w : 0.0;      /* found / limit */
            ws += sv; ws2 += sv * sv;
            double tt = w * secs; wsec += tt; wsec2 += tt * tt;
            if (k == knuth_mode || (k & (k - 1)) == 0) {
                double n = k, m = sw / n, se = sqrt(fmax(sw2 / n - m * m, 0) / n);
                double mz = wz / n, sez = sqrt(fmax(wz2 / n - mz * mz, 0) / n);
                double mt = wsec / n, set = sqrt(fmax(wsec2 / n - mt * mt, 0) / n);
                printf("walks %d leaves_est %.4g (se %.2g) reached %lld", k, m, se, nleaf);
                if (knuth_solve) printf(" | z4surv_est %.4g (se %.2g) found/limit_est %.3g | proc_secs_est %.4g (se %.2g)", mz, sez, ws / n, mt, set);
                printf(" | depth:");
                for (int L = 1; L < 7; L++) printf(" %.3g", wdepth[L] / n);
                if (knuth_pz4) {
                    double mp = wpass / n, sep = sqrt(fmax(wpass2 / n - mp * mp, 0) / n);
                    printf(" | pz4 kills:");
                    for (int L = 2; L < 7; L++) printf(" %.3g", kn_pzkill[L] / n);
                    printf(" secs:");
                    for (int L = 2; L < 7; L++) printf(" %.3g", kn_pzsec[L] / n);
                    printf(" leaves passing all pz4 %.4g (se %.2g)", mp, sep);
                }
                printf(" (%lds)%c", (long)(time(NULL) - t0), 10);
                fflush(stdout);
            }
        }
        return 0;
    }
    t_start = t_last = time(NULL);
    rows[0] = 0; nrows = 1;
    toponly = 1; dfs(parts, ncl); toponly = 0;              /* count top-level branches first */
    printf("top-level branches: %lld%c", top_total, 10); fflush(stdout);
    dfs(parts, ncl);
    snprintf(fn, sizeof fn, "%s.stats", argv[2]);
    FILE *fs = fopen(fn, "w");
    fprintf(fs, "{\"z4dead\": %lld, \"dup\": %lld, \"leaves\": %lld, \"processed\": %lld, \"empty\": %lld, \"root\": %lld, \"deeper\": %lld, \"heavy\": %lld, \"found\": %lld, \"limit\": %lld, \"m1_prunes\": %lld, \"nodes\": %lld, \"sample\": %g, \"part\": [%d, %d], \"seconds\": %ld}\n",
            st_z4, st_dup, st_leaves, st_proc, st_empty, st_root, st_deeper, st_heavy, st_found, st_limit, m1_prunes, st_nodes, samp, part_r, part_m, (long)(time(NULL) - t_start));
    fclose(fs); fclose(fsurv);
    if (fdbg) fclose(fdbg);
    printf("DONE leaves %lld processed %lld empty %lld root %lld deeper %lld heavy %lld found %lld limit %lld (%lds)\n",
           st_leaves, st_proc, st_empty, st_root, st_deeper, st_heavy, st_found, st_limit, (long)(time(NULL) - t_start));
    return 0;
}
