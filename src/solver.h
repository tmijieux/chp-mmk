#ifndef CHP_SOLVER_H
#define CHP_SOLVER_H

typedef enum chp_solver_type_ chp_solver_type;
typedef struct chp_solver_ chp_solver;

#include "equation.h"

#define SOLVER_EPSILON   1e-8
#define MAX_NB_ITER      20000

struct chp_solver_ {
    int Nx, Ny;
    double B, Cx, Cy;
    int (*solve)(int Nx, int Ny, const double B,
                 const double Cx, const double Cy,
                 double const *rhs, double *X, double *tmp[]);
    double *tmp[4];
};

void gemv_5diag(
    int const Nx, int const Ny, double const B,
    double const Cx, double const Cy,  double const *X, double *AX);
void chp_solver_init(chp_solver *S, chp_equation *eq, const char *solver_arg);
int chp_solver_run(chp_solver *S, double const *rhs, double const *X0, double *X);
void chp_solver_free(chp_solver *S);

#endif // CHP_SOLVER_H
