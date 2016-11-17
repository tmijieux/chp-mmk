#ifndef CHP_EQUATION_H
#define CHP_EQUATION_H

struct chp_equation;

#include "func.h"
#include "proc.h"

struct chp_equation {
    double Lx_min, Lx_max;
    double Ly_min, Ly_max;
    double dx, dy;
    double D, B, Cx, Cy;
    int Nx, Ny;

    double *X, *Y;
    double *top, *bottom, *left, *right;
    double *rhs;

    double *U0, *U1;
};

void chp_equation_grid_init(struct chp_equation *eq);
void chp_equation_init(struct chp_equation *eq,
                       int rank, int group_size, int r, int s);
void chp_equation_border_init(struct chp_proc *proc,
                              struct chp_equation *eq, struct chp_func *func);
void chp_equation_alloc(struct chp_equation *eq);


#define EQUATION_SPLAT_ARGS(eq) (eq)->Nx,(eq)->Ny,(eq)->B,(eq)->Cx,(eq)->Cy,(eq)->rhs,(eq)->U0,(eq)->U1

#endif // CHP_EQUATION_H
