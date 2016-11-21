#include <float.h>
#include "equation.h"
#include "util.h"
#include "proc.h"

void chp_equation_grid_init(struct chp_equation *eq)
{
    double dx = eq->dx;
    double dy = eq->dy;
    double *X = eq->X;
    double *Y = eq->Y;

    double m = DBL_MAX, M = DBL_MAX;
    for (int i = 0; i < eq->Nx; ++i) {
        X[i] = eq->Lx_min + (i+1)*dx;
        double l = fabs(X[i] - eq->prev_border_x);
        if (l < m) {
            m = l;
            eq->prev_border_col = i;
        }
        l = fabs(X[i] - eq->next_border_x);
        if (l < M) {
            M = l;
            eq->next_border_col = i;
        }
    }

    for (int i = 0; i < eq->Ny; ++i)
        Y[i] = eq->Ly_min + (i+1)*dy;
}

void chp_equation_init(
    struct chp_equation *eq, int rank, int group_size,
    int recouvr, int NNX, int NNY)
{
    int Nx = NNX / group_size;
    int Ny = NNY;
    if (rank < (NNX%group_size))
        ++ Nx;

    double recouvrD = ((double)recouvr/NNX);
    double Lx_min = max(0.0, ((double)rank / group_size) - recouvrD/2.0);
    double Lx_max = min(1.0, ((double)(rank+1) / group_size) + recouvrD/2.0);
    double Ly_min = 0.0, Ly_max = 1.0;

    if (rank > 0)
        eq->prev_border_x = Lx_min + recouvrD;
    if (rank < group_size-1)
        eq->next_border_x = Lx_max - recouvrD;

    //printf("r=%d :: Nx=%d :: [%g, %g]\n", rank, Nx, Lx_min, Lx_max);
    double dx = (Lx_max-Lx_min) / (Nx+1);
    double dy = 1.0 / (Ny+1);

    double D = 1.0;
    double B = 2 * D / SQUARE(dx) + 2 * D / SQUARE(dy);
    double Cx = -D / SQUARE(dx);
    double Cy = -D / SQUARE(dy);

    eq->Nx = Nx; eq->Ny = Ny;
    eq->Lx_min = Lx_min; eq->Lx_max = Lx_max;
    eq->Ly_min = Ly_min; eq->Ly_max = Ly_max;
    eq->dx = dx; eq->dy = dy;
    eq->D = D; eq->B = B; eq->Cx = Cx; eq->Cy = Cy;
}

void chp_equation_alloc(struct chp_equation *eq)
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
    eq->U0 = tdp_vector_new(N);
    eq->U1 = tdp_vector_new(N);
}

void chp_equation_free(struct chp_equation *eq)
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
    struct chp_proc *proc, struct chp_equation *eq, struct chp_func *func)
{
    (void) proc;

    int Nx = eq->Nx, Ny = eq->Ny;
    double *X = eq->X, *Y = eq->Y;

    if (func->type == CHP_STATIONARY)
        func->rhs(Nx, Ny, X, Y, eq->rhs);
    else
        func->rhs_unsta(Nx, Ny, X, Y, eq->rhs, 1.0, 1.0, 0.0);

    func->bottom(Nx, X, eq->bottom, eq->Ly_min);
    func->top(Nx, X, eq->top, eq->Ly_max);
    func->left(Ny, Y, eq->left, eq->Lx_min);
    func->right(Ny, Y, eq->right, eq->Lx_max);
    //func->U(Nx, Ny, X, Y, Lx_min, Lx_max, Ly, Uexact);
}


