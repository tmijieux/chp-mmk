#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "util.h"
#include "solver.h"
#include "equation.h"

static inline void gemv_5diag(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict X, double * restrict AX)
{
    {
        int const k = 0;
        AX[k] = B*X[k] + Cx*X[k+1] + Cy*X[k+Nx];
    }
    for (int i = 1; i < Nx-1; ++i)
    {
        int const k = i;
        AX[k] = B*X[k] + Cx*(X[k-1]+X[k+1]) + Cy*X[k+Nx];
    }
    {
        int const k = Nx-1;
        AX[k] = B*X[k] + Cx*X[k-1] + Cy*X[k+Nx];
    }

    // i = numéro de bloc-ligne de taille Nx | vertical
    for (int i = 1; i < Ny-1; ++i)
    {
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
    for (int i = 1; i < Nx-1; ++i)
    {
        int const k = Nx*Ny-Nx+i;
        AX[k] = B*X[k] + Cx*(X[k-1]+X[k+1]) + Cy*X[k-Nx];
    }
    {
        int const k = Nx*Ny-1;
        AX[k] = B*X[k] + Cx*X[k-1] + Cy*X[k-Nx];
    }
}

static int jacobi_5diag(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double *restrict X, double * restrict tmp[])
{
    int const N = Nx*Ny;
    double *Ax = tmp[0];
    double *X_old = X;
    double *X_new = tmp[1];

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        SWAP_POINTER(X_old, X_new);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int const k = j*Nx+i;
                const double b = rhs[k];
                const double b0 = ((j != 0)    ? (Cy * X_old[k-Nx]) : 0.0);
                const double b1 = ((i != 0)    ? (Cx * X_old[k-1])  : 0.0);
                const double b2 = ((i != Nx-1) ? (Cx * X_old[k+1])  : 0.0);
                const double b3 = ((j != Ny-1) ? (Cy * X_old[k+Nx]) : 0.0);
                X_new[k] = (b-(b0+b1+b2+b3)) / B;
            }
        }

        gemv_5diag(Nx, Ny, B, Cx, Cy, X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(SOLVER_EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

static inline void jacobi_5diag2_impl(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double *restrict X_new, double const *restrict X_old)
{
    // bottom
    {
        int const k = 0;
        X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*X_old[k+Nx])) / B;
    }
    for (int i = 1; i < Nx-1; ++i)
    {
        int const k = i;
        X_new[k] = (rhs[k] - (Cx*(X_old[k-1]+X_old[k+1]) + Cy*X_old[k+Nx])) / B;
    }
    {
        int const k = Nx-1;
        X_new[k] = (rhs[k] - (Cx*X_old[k-1] + Cy*X_old[k+Nx])) / B;
    }

    // middle
    for (int i = 1; i < Ny-1; ++i)
    {
        {
            int const k = Nx*i;
            X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*(X_old[k-Nx]+X_old[k+Nx]))) / B;
        }
        for (int j = 1; j < Nx-1; ++j)
        {
            int const k = Nx*i + j;
            X_new[k] = (rhs[k] - (Cx*(X_old[k-1]+X_old[k+1])
                                  + Cy*(X_old[k-Nx]+X_old[k+Nx]))) / B;
        }
        {
            int const k = Nx*i + Nx - 1;
            X_new[k] = (rhs[k] - (Cx*X_old[k-1] + Cy*(X_old[k-Nx]+X_old[k+Nx]))) / B;
        }
    }

    // top
    {
        int const k = Nx*Ny-Nx;
        X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*X_old[k-Nx])) / B;
    }
    for (int i = 1; i < Nx-1; ++i)
    {
        int const k = Nx*Ny-Nx+i;
        X_new[k] = (rhs[k] - (Cx*(X_old[k-1]+X_old[k+1]) + Cy*X_old[k-Nx])) / B;
    }
    {
        int const k = Nx*Ny-1;
        X_new[k] = (rhs[k] - (Cx*X_old[k-1] + Cy*X_old[k-Nx])) / B;
    }
}

static int jacobi_5diag2(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double * restrict X, double *restrict tmp[])
{
    int const N = Nx*Ny;
    double *Ax = tmp[0];
    double *X_old = X;
    double *X_new = tmp[1];

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        SWAP_POINTER(X_old, X_new);
        jacobi_5diag2_impl(Nx, Ny, B, Cx, Cy, rhs, X_new, X_old);
        gemv_5diag(Nx, Ny, B, Cx, Cy, X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(SOLVER_EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

static int gauss_seidel_5diag(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double * restrict X,  double *restrict tmp[])
{
    int const N = Nx*Ny;
    double *Ax = tmp[0];
    double *X_old = X;
    double *X_new = tmp[1];

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        SWAP_POINTER(X_old, X_new);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int const k = j*Nx+i;
                const double b = rhs[k];
                const double b0 = ((j != 0)    ? (Cy * X_new[k-Nx]) : 0.0);
                const double b1 = ((i != 0)    ? (Cx * X_new[k-1])  : 0.0);
                const double b2 = ((i != Nx-1) ? (Cx * X_old[k+1])  : 0.0);
                const double b3 = ((j != Ny-1) ? (Cy * X_old[k+Nx]) : 0.0);
                X_new[k] = (b-(b0+b1+b2+b3)) / B;
            }
        }

        gemv_5diag(Nx, Ny, B, Cx, Cy, X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(SOLVER_EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

static inline void gauss_seidel_5diag2_impl(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double *restrict X_new, double const *restrict X_old)
{
    // bottom
    {
        int const k = 0;
        X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*X_old[k+Nx])) / B;
    }
    for (int i = 1; i < Nx-1; ++i) {
        int const k = i;
        X_new[k] = (rhs[k] - (Cx*(X_new[k-1]+X_old[k+1]) + Cy*X_old[k+Nx])) / B;
    }
    {
        int const k = Nx-1;
        X_new[k] = (rhs[k] - (Cx*X_new[k-1] + Cy*X_old[k+Nx])) / B;
    }

    // middle
    for (int i = 1; i < Ny-1; ++i) {
        {
            int const k = Nx*i;
            X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*(X_new[k-Nx]+X_old[k+Nx]))) / B;
        }
        for (int j = 1; j < Nx-1; ++j) { // j = ligne au sein du bloc-ligne
            int const k = Nx*i + j;
            X_new[k] = (rhs[k] - (Cx*(X_new[k-1]+X_old[k+1])
                                  + Cy*(X_new[k-Nx]+X_old[k+Nx]))) / B;
        }
        {
            int const k = Nx*i + Nx - 1;
            X_new[k] = (rhs[k] - (Cx*X_new[k-1] + Cy*(X_new[k-Nx]+X_old[k+Nx]))) / B;
        }
    }

    // top
    {
        int const k = Nx*Ny-Nx;
        X_new[k] = (rhs[k] - (Cx*X_old[k+1] + Cy*X_new[k-Nx])) / B;
    }
    for (int i = 1; i < Nx-1; ++i) {
        int const k = Nx*Ny-Nx+i;
        X_new[k] = (rhs[k] - (Cx*(X_new[k-1]+X_old[k+1]) + Cy*X_new[k-Nx])) / B;
    }
    {
        int const k = Nx*Ny-1;
        X_new[k] = (rhs[k] - (Cx*X_new[k-1] + Cy*X_new[k-Nx])) / B;
    }
}

static int gauss_seidel_5diag2(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double * restrict X,  double * restrict tmp[])
{
    int const N = Nx*Ny;
    double *Ax = tmp[0];
    double *X_old = X;
    double *X_new = tmp[1];

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        SWAP_POINTER(X_old, X_new);
        gauss_seidel_5diag2_impl(Nx, Ny, B, Cx, Cy, rhs, X_new, X_old);
        gemv_5diag(Nx, Ny, B, Cx, Cy, X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(SOLVER_EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

static int CG_5diag(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const * restrict rhs, double * restrict X, double * restrict tmp[])
{
    int const N = Nx*Ny;
    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    double *Ax = tmp[0];
    double *P = tmp[1];
    double *Q = tmp[2];
    double *R = tmp[3];

    gemv_5diag(Nx, Ny, B, Cx, Cy, X, Ax);
    cblas_dcopy(N, rhs, 1, R, 1);
    cblas_daxpy(N, -1.0, Ax, 1, R, 1);
    cblas_dcopy(N, R, 1, P, 1);

    double gammaNew = cblas_ddot(N, R, 1, R, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        gemv_5diag(Nx, Ny, B, Cx, Cy, P, Q);
        double alpha = gammaNew / cblas_ddot(N, P, 1, Q, 1);

        cblas_daxpy(N, alpha, P, 1, X, 1);
        cblas_daxpy(N, -alpha, Q, 1, R, 1);

        double gammaOld = gammaNew;
        gammaNew = cblas_ddot(N, R, 1, R, 1);

        if ((gammaNew/nB) <= SQUARE(SOLVER_EPSILON))
            break;

        double beta = gammaNew / gammaOld;
        cblas_dscal(N, beta, P, 1);
        cblas_daxpy(N, 1.0, R, 1, P, 1);
    }

    return iter;
}

void chp_solver_init(chp_solver *S, chp_equation *eq, const char *solver_arg)
{
    memset(S, 0, sizeof*S);

    if (!strcmp(solver_arg, "CG"))
        S->solve = &CG_5diag;
    else if (!strcmp(solver_arg, "GS"))
        S->solve = &gauss_seidel_5diag2;
    else if (!strcmp(solver_arg, "J"))
        S->solve = &jacobi_5diag2;
    else {
        fprintf(stderr, "Solveur invalide: '%s'\n", solver_arg);
        exit(EXIT_FAILURE);
    }

    int N = eq->Nx * eq->Ny;
    for (int i = 0; i < 4; ++i)
        S->tmp[i] = tdp_vector_new(N);

    S->Nx = eq->Nx;
    S->Ny = eq->Ny;
    S->Cx = eq->Cx;
    S->Cy = eq->Cy;
    S->B = eq->B;
}

int chp_solver_run(chp_solver *S, double const *restrict rhs,
                   const double *X0, double *X)
{
    const int N = S->Nx * S->Ny;
    cblas_dcopy(N, X0, 1, X, 1);
    return S->solve(S->Nx, S->Ny, S->B, S->Cx, S->Cy, rhs, X, S->tmp);
}

void chp_solver_free(chp_solver *S)
{
    for (int i = 0; i < 4; ++i)
        free(S->tmp[i]);
    memset(S, 0, sizeof*S);
}
