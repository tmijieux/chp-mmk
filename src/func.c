#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "func.h"
#include "solver.h"

unsigned func_list_length = 0;
chp_func *func_list = NULL;

void
zero(int const N, double const *A, double *B, double v)
{
    (void) A;
    (void) v;
    memset(B, 0, N*sizeof B[0]);
}

static void
one(int const N, double const *A, double *B, double v)
{
    (void) A;
    (void) v;
    for (int i = 0; i < N; ++i)
        B[i] = 1.0;
}

static void
rhs_1(int const Nx, int const Ny, double const *X, double const *Y, double *F)
{
    for (int j = 0; j < Ny; ++j) {
        double yy = Y[j] - SQUARE(Y[j]);
        for (int i = 0; i < Nx; ++i)
            F[Nx*j+i] = 2*(yy + X[i] - SQUARE(X[i]));
    }
}

static void
U_1(int const Nx, int const Ny, double const *X, double const *Y,
    double const Lx, double const Ly, double *U)
{
    for (int j = 0; j < Ny; ++j) {
        double yy = Y[j] * (Ly-Y[j]);
        for (int i = 0; i < Nx; ++i)
            U[Nx*j+i] = X[i] * (Lx-X[i]) * yy;
    }
}

/*bottom OR top*/
static void
bottom_top_2(int const Nx, double const *X, double *G, double value)
{
    double const c1 = cos(value);
    for (int i = 0; i < Nx; ++i)
        G[i] = sin(X[i]) + c1;
}

/*right OR left*/
static void
right_left_2(int const Ny, double const *Y, double *G, double value)
{
    double const s1 = sin(value);
    for (int i = 0; i < Ny; ++i)
        G[i] = cos(Y[i]) + s1;
}

static void
rhs_2(int const Nx, int const Ny, double const *X, double const *Y, double *F)
{
    for (int j = 0; j < Ny; ++j) {
        double const s = cos(Y[j]);
        for (int i = 0; i < Nx; ++i)
            F[Nx*j+i] = sin(X[i]) + s;
    }
}

static void
U_2(int const Nx, int const Ny, double const *X, double const *Y,
    double const Lx, double const Ly, double *U)
{
    (void) Lx;
    (void) Ly;

    for (int j = 0; j < Ny; ++j) {
        double c = cos(Y[j]);
        for (int i = 0; i < Nx; ++i)
            U[Nx*j+i] = sin(X[i]) + c;
    }
}

static void
rhs_3(int const Nx, int const Ny, double const *X, double const *Y,
      double *F, double Lx, double Ly, double t)
{
    double c = cos(M_PI*t/2.0);

    for (int j = 0; j < Ny; ++j) {
        double e = exp(-SQUARE(Y[j] - 0.5*Ly)) * c;
        for (int i = 0; i < Nx; ++i)
            F[Nx*j+i] = exp(-SQUARE(X[i] - 0.5*Lx)) * e;
    }
}

REGISTER_FUNCTION(f=e^(-(x-(Lx/2)^2))*e^(-(y-(Ly/2)^2))*cos(PI/2*t),
                  UNSTATIONARY,
                  zero, zero,
                  one, one,
                  NULL, rhs_3, NULL);

REGISTER_FUNCTION(f=sin(x)+cos(y),
                  STATIONARY,
                  bottom_top_2, bottom_top_2,
                  right_left_2, right_left_2,
                  rhs_2, NULL, U_2);

REGISTER_FUNCTION(f=2*(x-x^2+y-y^2),
                  STATIONARY,
                  zero, zero,
                  zero, zero,
                  rhs_1, NULL, U_1);


static void get_func_by_id(unsigned idx, chp_func *func)
{
    if (idx >= func_list_length) {
        fprintf(stderr, "No such func!\n");
        exit(EXIT_FAILURE);

    }

    chp_func *l = func_list;
    while (idx--)
        l = l->next;

    *func = *l;
}

void chp_func_init(chp_func *func, unsigned idx, chp_proc *P)
{
    memset(func, 0, sizeof*func);
    get_func_by_id(idx, func);

    if (P->rank > 0)
        func->left = zero;

    if (P->rank < P->group_size-1)
        func->right = zero;
}

void chp_print_func_list(void)
{
    printf("Function list:\n");
    chp_func *list = func_list;
    int i = 0;
    while (list != NULL) {
        printf("%d: %s\n", i, list->name);
        list = list->next;
        ++ i;
    }

}
