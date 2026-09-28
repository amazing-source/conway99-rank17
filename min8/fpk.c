/* fpk.c: Fincke-Pohst for a positive definite integer Gram G (n <= 24): all x != 0 (up to sign) with x^T G x == target.
   Input (text): n target, then G row-major.  Output (text): one vector per line.  LLL is not applied (caller supplies a
   reasonably reduced basis); exact integer check at the leaves. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int n, target, G[24][24], x[24];
static double Qd[24][24], L[24][24];
static long long cnt = 0;
static void rec(int i, double rem, double *c) {
    /* x_i range given x_{i+1..n-1}: q_ii (x_i + sum_j q_ij x_j)^2 <= rem  (upper-triangular Cholesky form) */
    double ci = 0; for (int j = i + 1; j < n; j++) ci += Qd[i][j] * x[j];
    double r = sqrt(fmax(rem, 0) / Qd[i][i]);
    int lo = (int)ceil(-ci - r - 1e-9), hi = (int)floor(-ci + r + 1e-9);
    for (int v = lo; v <= hi; v++) {
        x[i] = v;
        double t = Qd[i][i] * (v + ci) * (v + ci);
        if (t > rem + 1e-9) continue;
        if (i == 0) {
            int f = -1; for (int k = 0; k < n; k++) if (x[k]) { f = k; break; }
            if (f < 0 || x[f] < 0) continue;
            long long s = 0; for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) s += (long long)x[a] * G[a][b] * x[b];
            if (s == target) { for (int k = 0; k < n; k++) printf("%d ", x[k]); printf("\n"); cnt++; }
        } else rec(i - 1, rem - t, c);
    }
    x[i] = 0;
}
int main(void) {
    if (scanf("%d %d", &n, &target) != 2) return 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (scanf("%d", &G[i][j]) != 1) return 1;
    /* q_ij from Cholesky: G = R^T R, q_ii = R_ii^2, q_ij = R_ij / R_ii */
    double R[24][24] = {{0}};
    for (int j = 0; j < n; j++) {
        double s = G[j][j]; for (int k = 0; k < j; k++) s -= R[k][j] * R[k][j];
        R[j][j] = sqrt(s);
        for (int i = j + 1; i < n; i++) { double t = G[j][i]; for (int k = 0; k < j; k++) t -= R[k][j] * R[k][i]; R[j][i] = t / R[j][j]; }
    }
    for (int i = 0; i < n; i++) { Qd[i][i] = R[i][i] * R[i][i]; for (int j = i + 1; j < n; j++) Qd[i][j] = R[i][j] / R[i][i]; }
    for (int i = 0; i < n; i++) x[i] = 0;
    rec(n - 1, target, NULL);
    fprintf(stderr, "count %lld\n", cnt);
    return 0;
}
