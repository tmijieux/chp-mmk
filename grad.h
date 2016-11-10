#ifndef PROJCHP_GRAD_H
#define PROJCHP_GRAD_H

void matrix_5diag_conjugate_gradient(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *rhs, double *X);

void matrix_5diag_sym_product(
    int const Nx, int const Ny,
    double const B, double const Cx, double const Cy,
    double const *X, double *AX);

void vector_compute_RHS(int const Nx, int const Ny,
                        double const Cx, double const Cy,
                        double const *h, double const *g,
                        double */*inout*/RHS);

#endif // PROJCHP_GRAD_H
