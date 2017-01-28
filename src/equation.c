#include <float.h>
#include <mpi.h>

#include "equation.h"
#include "util.h"
#include "proc.h"
#include "cblas.h"

static void init_grid(chp_equation *eq)
{
    double dx = eq->dx;
    double dy = eq->dy;
    double *X = eq->X;
    double *Y = eq->Y;

    double m1 = DBL_MAX, M1 = DBL_MAX;
    double m2 = DBL_MAX, M2 = DBL_MAX;
    for (int i = 0; i < eq->Nx; ++i) {
        X[i] = eq->Lx_min + (i+1)*dx;
        double l = fabs(X[i] - eq->prev_border_x);
        if (l < m1) {
            m2 = m1;
            eq->prev_border_col2 = eq->prev_border_col;

            m1 = l;
            eq->prev_border_col = i;
        } else if (l < m2) {
            m2 = l;
            eq->prev_border_col2 = i;
        }

        l = fabs(X[i] - eq->next_border_x);
        if (l < M1) {
            M2 = M1;
            eq->next_border_col2 = eq->next_border_col;
            M1 = l;
            eq->next_border_col = i;
        } else if (l < M2) {
            M2 = l;
            eq->next_border_col2 = i;
        }
    }
    if (eq->next_border_col2 < eq->next_border_col)
        SWAP_VARS(eq->next_border_col2, eq->next_border_col);

    if (eq->prev_border_col2 < eq->prev_border_col)
        SWAP_VARS(eq->prev_border_col2, eq->prev_border_col);

    for (int i = 0; i < eq->Ny; ++i)
        Y[i] = eq->Ly_min + (i+1)*dy;
}

static void init_prop(chp_equation *eq, int rank, int group_size,
                      int NNX, int NNY, int recouvr,
                      double Lx, double Ly, double Tmax, int Nit, bool stationary)
{
    int Nx = NNX / group_size;
    int Ny = NNY;

        if (rank < (NNX % group_size))
        ++ Nx;

    double recouvrD = ((double)recouvr/NNX)*Lx;
    double Lx_min = max(0.0, ((double)rank / group_size)*Lx - recouvrD/2.0);
    double Lx_max = min(Lx, ((double)(rank+1) / group_size)*Lx + recouvrD/2.0);
    double Ly_min = 0.0, Ly_max = Ly;

    if (rank > 0)
        eq->prev_border_x = Lx_min + recouvrD;
    if (rank < group_size-1)
        eq->next_border_x = Lx_max - recouvrD;

    //printf("r=%d :: Nx=%d :: [%g, %g]\n", rank, Nx, Lx_min, Lx_max);
    double dx = (Lx_max-Lx_min) / (Nx+1);
    double dy = (Ly_max-Ly_min) / (Ny+1);

    double D = 1.0;
    double B = 2 * D / SQUARE(dx) + 2 * D / SQUARE(dy);
    double Cx = -D / SQUARE(dx);
    double Cy = -D / SQUARE(dy);

    eq->Nx = Nx; eq->Ny = Ny;
    eq->N = Nx*Ny;
    eq->Lx_min = Lx_min; eq->Lx_max = Lx_max;
    eq->Ly_min = Ly_min; eq->Ly_max = Ly_max;
    eq->dx = dx; eq->dy = dy;
    eq->D = D; eq->B = B; eq->Cx = Cx; eq->Cy = Cy;
    eq->Lx = Lx; eq->Ly = Ly;
    eq->Tmax = Tmax;
    eq->Nit = Nit;

    if (!stationary) {
        eq->dt = Tmax / Nit;
        eq->B += 1.0 / eq->dt;
    }

}

static void equation_alloc(chp_equation *eq)
{
    int Nx = eq->Nx, Ny = eq->Ny;
    int N = Nx*Ny;

    eq->top = tdp_vector_new(Nx);
    eq->bottom = tdp_vector_new(Nx);
    eq->right = tdp_vector_new(Ny);
    eq->left = tdp_vector_new(Ny);

    eq->X = tdp_vector_new(Nx);
    eq->Y = tdp_vector_new(Ny);
    //eq->Y = NULL;

    eq->rhs = tdp_vector_new(N);
    eq->rhs_f = tdp_vector_new(N);
    eq->U0 = tdp_vector_new(N);
    eq->U1 = tdp_vector_new(N);
}

void chp_equation_init(
    chp_equation *eq, chp_proc *P, struct gengetopt_args_info *opt, bool stationary)
{
    memset(eq, 0, sizeof*eq);
    init_prop(eq, P->rank, P->group_size,
              opt->resolutionX_arg, opt->resolutionY_arg,
              opt->recouvr_arg, opt->Lx_arg, opt->Ly_arg,
              opt->Tmax_arg, opt->Nit_arg, stationary);

    equation_alloc(eq);
    init_grid(eq);
}

void chp_equation_free(chp_equation *eq)
{
    free(eq->top);
    free(eq->bottom);
    free(eq->right);
    free(eq->left);
    free(eq->X);
    free(eq->Y);
    free(eq->rhs);
    free(eq->U0);
    free(eq->U1);
}

void chp_equation_border_init(
    chp_equation *eq, chp_func *func)
{

    int Nx = eq->Nx, Ny = eq->Ny;
    double *X = eq->X, *Y = eq->Y;

    func->bottom(Nx, X, eq->bottom, eq->Ly_min);
    func->top(Nx, X, eq->top, eq->Ly_max);

    func->left(Ny, Y, eq->left, eq->Lx_min);
    func->right(Ny, Y, eq->right, eq->Lx_max);
    //func->U(Nx, Ny, X, Y, Lx_min, Lx_max, Ly, Uexact);
}

void chp_equation_rhs_init(
    chp_equation *eq, chp_func *func, double t)
{
    if (func->type == CHP_STATIONARY)
        func->rhs(eq->Nx, eq->Ny, eq->X, eq->Y, eq->rhs_f);
    else {
        func->rhs_unsta(eq->Nx, eq->Ny, eq->X, eq->Y,
                        eq->rhs_f, eq->Lx, eq->Ly, t);
        cblas_daxpy(eq->N, 1.0/eq->dt, eq->U0, 1, eq->rhs_f, 1);
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
void chp_equation_apply_border_cond_RHS(chp_equation *eq)
{
    int Nx = eq->Nx, Ny = eq->Ny;

    cblas_dcopy(eq->N, eq->rhs_f, 1, eq->rhs, 1);

    cblas_daxpy(Nx, -eq->Cy, eq->bottom, 1, eq->rhs, 1);
    cblas_daxpy(Nx, -eq->Cy, eq->top, 1, eq->rhs+Nx*(Ny-1), 1);

    cblas_daxpy(Ny, -eq->Cx, eq->left, 1, eq->rhs, Nx);
    cblas_daxpy(Ny, -eq->Cx, eq->right, 1, eq->rhs+Nx-1, Nx);
}
