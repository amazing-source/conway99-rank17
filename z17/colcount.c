/* colcount.c — count the admissible COLUMNS of the w-matrix for V = Z^17 (structures, 28/09).
   Column j of the 42 x 17 matrix of the w's: entries w_{cd,s,j}, (c,d) the 21 cells, s = +/-.
   - odd exactly on the cut edges of a bipartition P | P^c of the 7 vertices (|P| = p in {1,2,3}): (w+,w-) in {+-1}^2;
   - on the other cells (w+,w-) in {0,+-2}^2;
   - squared norm sum = T (T = 28 (P_K)_jj: 28 for y-free coordinates, 24 or 20 otherwise);
   - stars: at each vertex v, sum_{d>v} a_vd + sum_{d<v} b_dv = 0 with a = (w+ + w-)/2, b = (w+ - w-)/2.
   Prints the number of admissible columns per (p, T). Usage: colcount */
#include <stdio.h>
#include <string.h>
static int CC[21], CD[21], cut[21], T, star[7];
static long long cnt;
static int lastcell[7];
static void rec(int k, int norm) {
    if (norm > T) return;
    if (k == 21) { if (norm == T) cnt++; return; }
    int c = CC[k], d = CD[k];
    static const int oddab[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    static const int evab[9][3] = {{0, 0, 0}, {1, 1, 4}, {-1, -1, 4}, {1, -1, 4}, {-1, 1, 4}, {2, 0, 8}, {-2, 0, 8}, {0, 2, 8}, {0, -2, 8}};
    int no = cut[k] ? 4 : 9;
    for (int o = 0; o < no; o++) {
        int a, b, dn;
        if (cut[k]) { a = oddab[o][0]; b = oddab[o][1]; dn = 2; } else { a = evab[o][0]; b = evab[o][1]; dn = evab[o][2]; }
        star[c] += a; star[d] += b;
        int ok = 1;
        if (lastcell[c] == k && star[c]) ok = 0;
        if (lastcell[d] == k && star[d]) ok = 0;
        if (ok) rec(k + 1, norm + dn);
        star[c] -= a; star[d] -= b;
    }
}
int main(void) {
    int k = 0;
    for (int c = 0; c < 7; c++) for (int d = c + 1; d < 7; d++) { CC[k] = c; CD[k] = d; k++; }
    for (int v = 0; v < 7; v++) lastcell[v] = -1;
    for (int x = 0; x < 21; x++) { lastcell[CC[x]] = x; lastcell[CD[x]] = x; }
    int Ts[3] = {28, 24, 20};
    for (int p = 1; p <= 3; p++) {
        for (int x = 0; x < 21; x++) cut[x] = ((CC[x] < p) != (CD[x] < p));
        for (int ti = 0; ti < 3; ti++) {
            T = Ts[ti]; cnt = 0; memset(star, 0, sizeof star);
            rec(0, 0);
            printf("p=%d T=%d columns=%lld\n", p, T, cnt);
            fflush(stdout);
        }
    }
    return 0;
}
