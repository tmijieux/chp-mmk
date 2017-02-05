#include <cstdlib>
#include <cstring>

#include "cblas.h"
#include "util.h"
#include "solver.hpp"
#include "equation.hpp"
#include "error.hpp"


namespace chp {

void jacobi_solver::impl(const double * __restrict__ rhs,
                         double *__restrict__ X_new,
                         const double *__restrict__ X_old)
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

int jacobi_solver::_solve(const double *__restrict__ rhs, double *__restrict__ X)
{
    double *Ax = &_tmp[0][0];
    double *X_new = &_tmp[1][0];
    double *X_old = X;

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        std::swap(X_old, X_new);
        impl(rhs, X_new, X_old);
        gemv_5diag(X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

void solver::gemv_5diag(const double * __restrict__ X, double * __restrict__ AX)
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
void gauss_seidel_solver::impl(const double * __restrict__ rhs,
                               double *__restrict__ X_new,
                               const double *__restrict__ X_old)
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

int gauss_seidel_solver::_solve(const double * __restrict__ rhs, double * __restrict__ X)
{
    int const N = Nx*Ny;
    double *Ax = &_tmp[0][0];
    double *X_old = X;
    double *X_new = &_tmp[1][0];

    double nB = cblas_ddot(N, rhs, 1, rhs, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        std::swap(X_old, X_new);
        impl(rhs, X_new, X_old);
        gemv_5diag(X_new, Ax);
        cblas_daxpy(N, -1.0, rhs, 1, Ax, 1);
        double gamma = cblas_ddot(N, Ax, 1, Ax, 1);
        if ((gamma/nB) <= SQUARE(EPSILON))
            break;
    }

    if (X != X_new)
        cblas_dcopy(N, X_new, 1, X, 1);

    return iter;
}

int CG_solver::_solve(const double * __restrict__ rhs,  double * __restrict__ X)
{
    int const N = Nx*Ny;
    double nB = cblas_ddot(N, rhs, 1, rhs, 1);

    double *Ax = &_tmp[0][0];
    double *P = &_tmp[1][0];
    double *Q = &_tmp[2][0];
    double *R = &_tmp[3][0];

    gemv_5diag(X, Ax);
    cblas_dcopy(N, rhs, 1, R, 1);
    cblas_daxpy(N, -1.0, Ax, 1, R, 1);
    cblas_dcopy(N, R, 1, P, 1);

    double gammaNew = cblas_ddot(N, R, 1, R, 1);
    int iter;
    for (iter = 0; iter < MAX_NB_ITER; ++iter) {
        gemv_5diag(P, Q);
        double alpha = gammaNew / cblas_ddot(N, P, 1, Q, 1);

        cblas_daxpy(N, alpha, P, 1, X, 1);
        cblas_daxpy(N, -alpha, Q, 1, R, 1);

        double gammaOld = gammaNew;
        gammaNew = cblas_ddot(N, R, 1, R, 1);

        if ((gammaNew/nB) <= SQUARE(EPSILON))
            break;

        double beta = gammaNew / gammaOld;
        cblas_dscal(N, beta, P, 1);
        cblas_daxpy(N, 1.0, R, 1, P, 1);
    }
    return iter;
}

solver *solver::make_solver(equation const& eq, std::string const& solver_arg)
{

    if (solver_arg == "CG")
        return new CG_solver(eq);
    else if (solver_arg == "GS")
        return new gauss_seidel_solver(eq);
    else if (solver_arg == "J")
        return new jacobi_solver(eq);
    else
        throw exception("Solveur invalide: '%s'") % solver_arg;
}

solver::solver(equation const& eq):
    Nx(eq.Nx), Ny(eq.Ny), N(Nx*Ny), B(eq.B), Cx(eq.Cx), Cy(eq.Cy)
{
    for (int i = 0; i < 4; ++i)
        _tmp[i].resize(N);
}

int solver::solve(const double *__restrict__ rhs, const double *X0, double *X)
{
    cblas_dcopy(N, X0, 1, X, 1);
    return _solve(rhs, X);
}

}; // namespace chp
