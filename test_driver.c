#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

#include "error.h"
#include "perf/perf.h"

#include "grad.h"
#include "util.h"

static void
test_vector_compute_RHS(void)
{
    int const Nx = 10;
    int const Ny = 10;
    double Cx = 0.1, Cy = 0.1;

    struct chp_equation eq;
    eq.rhs = tdp_vector_new(Nx*Ny);
    eq.right = tdp_vector_new(Ny);
    eq.left = tdp_vector_new(Ny);
    eq.bottom = tdp_vector_new(Nx);
    eq.top = tdp_vector_new(Nx);

    tdp_vector_one(Nx, 3.0, eq.bottom);
    tdp_vector_one(Nx, 3.0, eq.top);
    tdp_vector_one(Ny, 1.0, eq.right);
    tdp_vector_one(Ny, 1.0, eq.left);
    
    eq.Nx = Nx; eq.Ny = Ny;
    eq.Cx = Cx; eq.Cy = Cy;
    

    vector_compute_RHS(&eq);
    tdp_vector_print(100, eq.rhs, stdout);
}

static void
test_matrix_5diag_sym_product(void)
{
    int const Nx = 10;
    int const Ny = 10;
    int const N = Nx*Ny;

    double *X = tdp_vector_new(N);
    double *AX = tdp_vector_new(N);
    double B = 3.0, Cx = 1.0, Cy = 1.0;

    for (int i = 0; i < N; ++i)
        X[i] = 1.0 * (i+1);

    matrix_5diag_sym_product(Nx, Ny, B, Cx, Cy, X, AX);
    tdp_vector_print(N, AX, stdout);
}

int main()
{
    test_vector_compute_RHS();
    test_matrix_5diag_sym_product();
    return EXIT_SUCCESS;
}

