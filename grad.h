#ifndef PROJCHP_GRAD_H
#define PROJCHP_GRAD_H

#include "equation.h"


#define EPSILON   1e-5
#define SWAP_POINTER(p1, p2)                     \
    do {                                        \
        void *tmp_MAXCRO__ = p1;                \
        p1 = p2;                                \
        p2 = tmp_MAXCRO__;                      \
    } while(0)                                  \


void matrix_5diag_jacobi(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X0, double *X);

void matrix_5diag_gauss_seidel(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X0, double *X);

void matrix_5diag_conjugate_gradient(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X);

void matrix_5diag_sym_product(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *X, double *AX);

void vector_compute_RHS(struct chp_equation *eq);

#endif // PROJCHP_GRAD_H
