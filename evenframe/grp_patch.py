p='evenframe.c'
s=open(p).read()
# ---- globals
s=s.replace('static int dynstart = 12;','''static int dynstart = 12;
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
}''')
# ---- canonicity test before the PSD work
old='''    const double *p = pr[x] + (size_t)t * n;
    double R2[MAXN][MAXN];'''
new='''    if (ngr > 1 && k < dynstart) {
        int sv = chosen[o]; chosen[o] = t; asg[o] = 1;
        int okc = canon_ok();
        asg[o] = 0; chosen[o] = sv;
        if (!okc) { canpruned++; return 0; }
    }
    const double *p = pr[x] + (size_t)t * n;
    double R2[MAXN][MAXN];'''
assert old in s; s=s.replace(old,new,1)
old='''    /* commit */
    chosen[o] = t;'''
new='''    /* commit */
    chosen[o] = t; asg[o] = 1;'''
assert old in s; s=s.replace(old,new,1)
old='''    /* undo */'''
new='''    /* undo */
    asg[o] = 0;'''
assert old in s; s=s.replace(old,new,1)
# ---- read generators of the labelled group, closure, maps
old='''    fclose(f);
    for (int g = 0; g < ngen; g++) {'''
new='''    int ngg = 0; short (*GG)[MAXN][MAXN] = NULL; unsigned char (*GGL)[7] = NULL;
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
    for (int g = 0; g < ngen; g++) {'''
assert old in s; s=s.replace(old,new,1)
# the budgets are read before the first fclose: move the group read after budgets -> the first replace anchored on the
# fclose following the budget section, which is what we want (order in file: gens, budgets, group)
s=s.replace('printf("%s %lld sym %lld bud %lld cn %lld |", r == 1','printf("%s %lld sym %lld bud %lld cn %lld can %lld grp %d |", r == 1')
s=s.replace('''"EXCLUDED"), nodes, sympruned, budpruned, cnpruned);''','''"EXCLUDED"), nodes, sympruned, budpruned, cnpruned, canpruned, ngr);''')
open(p,'w',newline='\n').write(s)
