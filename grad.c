#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "util.h"
#include "grad.h"

#define MAX_NB_ITER 20000
#define EPSILON   1e-5

void matrix_5diag_conjugate_gradient(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X)
{
    int const N = Nx*Ny;

    double *P = tdp_vector_new(N);
    double *Q = tdp_vector_new(N);
    double *R = tdp_vector_new(N);
    double *Ax = tdp_vector_new(N);

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
    cblas_dcopy(N, rhs, 1, R, 1);
    cblas_daxpy(N, -1.0, Ax, 1, R, 1);
    cblas_dcopy(N, R, 1, P, 1);
    free(Ax);

    double gammaNew = cblas_ddot(N, R, 1, R, 1);
    int k;
    for (k = 0; k < MAX_NB_ITER; ++k) {
        matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, P, Q);
        double alpha = gammaNew / cblas_ddot(N, P, 1, Q, 1);

        cblas_daxpy(N, alpha, P, 1, X, 1);
        cblas_daxpy(N, -alpha, Q, 1, R, 1);

        double gammaOld = gammaNew;
        gammaNew = cblas_ddot(N, R, 1, R, 1);

        if ((gammaNew / nB) <= SQUARE(EPSILON))
            break;

        double beta = gammaNew / gammaOld;
        cblas_dscal(N, beta, P, 1);
        cblas_daxpy(N, 1.0, R, 1, P, 1);
    }
    free(P);
    free(Q);
    free(R);
}

void matrix_5diag_sym_product(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *X, double *AX)
{
    {
        int const k = 0;
        AX[k] = B*X[k] + Cx*X[k+1] + Cy*X[k+Nx];
    }
    for (int i = 1; i < Nx-1; ++i) {
        int const k = i;
        AX[k] = B*X[k] + Cx*(X[k-1]+X[k+1]) + Cy*X[k+Nx];
    }
    {
        int const k = Nx-1;
        AX[k] = B*X[k] + Cx*X[k-1] + Cy*X[k+Nx];
    }

    // i = numéro de bloc-ligne de taille Nx | vertical
    for (int i = 1; i < Ny-1; ++i) {
        {
            int const k = Nx*i;
            AX[k] = B*X[k] + Cx*X[k+1] + Cy*(X[k-Nx] + X[k+Nx]);
        }
        for (int j = 1; j < Nx-1; ++j) { // j = ligne au sein du bloc-ligne
            int const k = Nx*i + j;
            AX[k] = B*X[k] + Cx*(X[k-1] + X[k+1]) + Cy*(X[k-Nx] + X[k+Nx]);
        }
        {
            int const k = Nx*i + Nx - 1;
            AX[k] = B*X[k] + Cx*X[k-1] + Cy*(X[k-Nx] + X[k+Nx]);
        }
    }

    {
        int const k = Nx*Ny-Nx;
        AX[k] = B*X[k] + Cx*X[k+1] + Cy*X[k-Nx];
    }
    for (int i = 1; i < Nx-1; ++i) {
        int const k = Nx*Ny-Nx+i;
        AX[k] = B*X[k] + Cx*(X[k-1]+X[k+1]) + Cy*X[k-Nx];
    }
    {
        int const k = Nx*Ny-1;
        AX[k] = B*X[k] + Cx*X[k-1] + Cy*X[k-Nx];
    }
}

/**
 * This method apply border condition on equation RHS
 *
 * g is condition on x border | size 2*Nx
 * (first half is bottom, second half is top)
 *
 * h is condition on y border | size 2*Ny
 * (first half is left, second half is right)
 *
 */
void vector_compute_RHS(int const Nx, int const Ny,
                        double const Cx, double const Cy,
                        double const *h, double const *g,
                        double */*inout*/RHS)
{
    cblas_daxpy(Nx, -Cy, g, 1, RHS, 1);
    cblas_daxpy(Ny, -Cx, h, 1, RHS, Nx);
    cblas_daxpy(Ny, -Cx, h+Ny, 1, RHS+Nx-1, Nx);
    cblas_daxpy(Nx, -Cy, g+Nx, 1, RHS+Ny*(Nx-1), 1);
}
