#include <float.h>
#include <mpi.h>

#include "equation.hpp"
#include "proc.hpp"

#include "util.h"
#include "cblas.h"

using namespace chp;

void equation::init_grid()
{
    double m1 = DBL_MAX, M1 = DBL_MAX;
    double m2 = DBL_MAX, M2 = DBL_MAX;
    for (int i = 0; i < Nx; ++i) {
        X[i] = Lx_min + (i+1)*dx;
        double l = fabs(X[i] - prev_border_x);
        if (l < m1) {
            m2 = m1;
            prev_border_col2 = prev_border_col;

            m1 = l;
            prev_border_col = i;
        } else if (l < m2) {
            m2 = l;
            prev_border_col2 = i;
        }

        l = fabs(X[i] - next_border_x);
        if (l < M1) {
            M2 = M1;
            next_border_col2 = next_border_col;
            M1 = l;
            next_border_col = i;
        } else if (l < M2) {
            M2 = l;
            next_border_col2 = i;
        }
    }
    if (next_border_col2 < next_border_col)
        std::swap(next_border_col2, next_border_col);

    if (prev_border_col2 < prev_border_col)
        std::swap(prev_border_col2, prev_border_col);

    for (int i = 0; i < Ny; ++i)
        Y[i] = Ly_min + (i+1)*dy;
}

void equation::init_prop(proc const& P, int NNX, int NNY, int recouvr, bool stationary)
{
    Nx = NNX / P.size();
    Ny = NNY;
    if (P.rank() < (NNX % P.size()))
        ++ Nx;
    N = Nx*Ny;

    double recouvrD = ((double)recouvr/NNX)*Lx;
    Lx_min = std::max(0.0, ((double)P.rank() / P.size())*Lx - recouvrD/2.0);
    Lx_max = std::min(Lx, ((double)(P.rank()+1) / P.size())*Lx + recouvrD/2.0);
    Ly_min = 0.0;
    Ly_max = Ly;

    if (P.rank() > 0)
        prev_border_x = Lx_min + recouvrD;
    if (P.rank() < P.size()-1)
        next_border_x = Lx_max - recouvrD;

    //printf("r=%d :: Nx=%d :: [%g, %g]\n", P.rank(), Nx, Lx_min, Lx_max);
    dx = (Lx_max-Lx_min) / (Nx+1);
    dy = (Ly_max-Ly_min) / (Ny+1);

    D = 1.0;
    B = 2 * D / SQUARE(dx) + 2 * D / SQUARE(dy);
    Cx = -D / SQUARE(dx);
    Cy = -D / SQUARE(dy);

    ASSERT_MSG( Nit > 0, "Le nombre de pas de temps doit être positif.");

    if (!stationary) {
        dt = Tmax / Nit;
        B += 1.0 / dt;
    }
}

void equation::alloc()
{
    top = tdp_vector_new(Nx);
    bottom = tdp_vector_new(Nx);
    right = tdp_vector_new(Ny);
    left = tdp_vector_new(Ny);

    X = tdp_vector_new(Nx);
    Y = tdp_vector_new(Ny);

    rhs = tdp_vector_new(N);
    rhs_f = tdp_vector_new(N);
    U0 = tdp_vector_new(N);
    U1 = tdp_vector_new(N);
}

equation::equation(proc const& P, struct gengetopt_args_info *opt, bool stationary):
    Lx(opt->Lx_arg), Ly(opt->Ly_arg), Nit(opt->Nit_arg), Tmax(opt->Tmax_arg)
{
    init_prop(P, opt->resolutionX_arg, opt->resolutionY_arg, opt->recouvr_arg, stationary);
    std::cerr << "DEBUG: " << Nx << std::endl;
    alloc();

    init_grid();
}

equation::~equation()
{
    free(top);
    free(bottom);
    free(right);
    free(left);
    free(X);
    free(Y);
    free(rhs);
    free(rhs_f);
    free(U0);
    free(U1);
}

void equation::border_init(func const &f)
{
    f._bottom(Nx, X, bottom, Ly_min);
    f._top(Nx, X, top, Ly_max);

    f._left(Ny, Y, left, Lx_min);
    f._right(Ny, Y, right, Lx_max);
    //f._U(Nx, Ny, X, Y, Lx_min, Lx_max, Ly, Uexact);
}

void equation::rhs_init(func const &f, double t)
{
    if (f._type == func::STATIONARY)
        f._rhs(Nx, Ny, X, Y, rhs_f);
    else {
        f._rhs_unsta(Nx, Ny, X, Y, rhs_f, Lx, Ly, t);
        cblas_daxpy(N, 1.0/dt, U0, 1, rhs_f, 1);
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
void equation::apply_border_cond_RHS()
{
    cblas_dcopy(N, rhs_f, 1, rhs, 1);

    cblas_daxpy(Nx, -Cy, bottom, 1, rhs, 1);
    cblas_daxpy(Nx, -Cy, top, 1, rhs+Nx*(Ny-1), 1);

    cblas_daxpy(Ny, -Cx, left, 1, rhs, Nx);
    cblas_daxpy(Ny, -Cx, right, 1, rhs+Nx-1, Nx);
}
