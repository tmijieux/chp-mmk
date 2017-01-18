#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "util.h"
#include "grad.h"
#include "equation.h"

#define MAX_NB_ITER 20000

#define PASTE2_(x,y) x##y
#define PASTE2(x,y) PASTE2_(x, y)


void matrix_5diag_jacobi(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X0, double *X)
{
    int const N = Nx*Ny;

    double *Ax = tdp_vector_new(N);
    cblas_dcopy(N, X0, 1, X, 1);
    matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
    cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);

    double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    for (int iter = 0; iter < MAX_NB_ITER; ++iter) {

        SWAP_POINTER(X0, X);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int const k = j*Nx + i;
                double b = rhs[k];
                b = b - ((j != 0) ? (Cy * X0[k-Nx]) : 0);
                b = b - ((i != 0) ? (Cx * X0[k-1]) : 0);
                b = b - ((i != Nx-1) ? (Cx * X0[k+1]) : 0);
                b = b - ((j != Ny-1) ? (Cy * X0[k+Nx]) : 0);
                X[k] = b / B;
            }
        }

        matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(EPSILON))
            break;
    }
    free(Ax);
}

void matrix_5diag_gauss_seidel(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X0, double *X)
{
    int const N = Nx*Ny;

    double *Ax = tdp_vector_new(N);
    cblas_dcopy(N, X0, 1, X, 1);
    matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
    cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);

    double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    for (int iter = 0; iter < MAX_NB_ITER; ++iter) {
      
        SWAP_POINTER(X0, X);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int const k = j*Nx+i;
                double b = rhs[k];
                b = b - ((j != 0) ? (Cy * X[k-Nx]) : 0.0);
                b = b - ((i != 0) ? (Cx * X[k-1]) : 0.0);
                b = b - ((i != Nx-1) ? (Cy * X0[k+1]) : 0.0);
                b = b - ((j != Ny-1) ? (Cy * X0[k+Nx]) : 0.0);
                X[k] = b / B;
            }
        }

	matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(EPSILON))
            break;
    }
    free(Ax);
}

/*void matrix_5diag_gauss_seidel(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X0, double *X)
{
    int const N = Nx*Ny;

    double *Ax = tdp_vector_new(N);
    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    for (int iter = 0; iter < MAX_NB_ITER; ++iter) {
        SWAP_POINTER(X0, X);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int const k = j*Nx+i;
                double b = rhs[k];
                b = b - ((j != 0) ? (Cy * X[k-Nx]) : 0.0);
                b = b - ((i != 0) ? (Cx * X[k-1]) : 0.0);
                b = b - ((i != Nx-1) ? (Cy * X0[k+1]) : 0.0);
                b = b - ((j != Ny-1) ? (Cy * X0[k+Nx]) : 0.0);
                X[k] = b / B;
            }
        }
        matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(EPSILON))
            break;
    }
    free(Ax);
}*/

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
void vector_compute_RHS(struct chp_equation *eq)
 {
    int Nx = eq->Nx, Ny = eq->Ny;

    cblas_daxpy(Nx, -eq->Cy, eq->bottom, 1, eq->rhs, 1);
    cblas_daxpy(Nx, -eq->Cy, eq->top, 1, eq->rhs+Nx*(Ny-1), 1);

    cblas_daxpy(Ny, -eq->Cx, eq->left, 1, eq->rhs, Nx);
    cblas_daxpy(Ny, -eq->Cx, eq->right, 1, eq->rhs+Nx-1, Nx);
}
