#ifndef CHP_EQUATION_H
#define CHP_EQUATION_H

#include <stdbool.h>

typedef struct chp_equation_ chp_equation;

#include "func.h"
#include "proc.h"
#include "cmdline.h"

struct chp_equation_ {
    double Lx_min, Lx_max;
    double Ly_min, Ly_max;
    double Lx, Ly;

    int Nit;
    double Tmax;
    double dx, dy, dt;
    double D, B, Cx, Cy;
    int Nx, Ny;
    int N;

    double *X, *Y;
    double *top, *bottom, *left, *right;
    double *rhs, *rhs_f;

    double *U0, *U1;

    double next_border_x, prev_border_x;
    int next_border_col, prev_border_col;
    int next_border_col2, prev_border_col2;
};


void chp_equation_init(chp_equation *eq, chp_proc *P,
                       struct gengetopt_args_info *opt, bool stationary);
void chp_equation_grid_init(chp_equation *eq);
void chp_equation_alloc(chp_equation *eq);
void chp_equation_border_init(chp_equation *eq, chp_func *func);
void chp_equation_rhs_init(chp_equation *eq, chp_func *func, double t);
void chp_equation_free(chp_equation *eq);
void chp_equation_apply_border_cond_RHS(chp_equation *eq);

#define EQUATION_SPLAT_ARGS(eq) \
    (eq)->Nx,(eq)->Ny,(eq)->B,(eq)->Cx,(eq)->Cy,(eq)->rhs,(eq)->U0,(eq)->U1

#endif // CHP_EQUATION_H
