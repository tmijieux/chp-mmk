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
    double *A = tdp_vector_new(100);
    int const Nx = 10;
    int const Ny = 10;
    double Cx = 0.1, Cy = 0.1;

    double *g = tdp_vector_new(2*Nx);
    double *h = tdp_vector_new(2*Ny);

    tdp_vector_one(2*Nx, 3.0, g);
    tdp_vector_one(2*Ny, 1.0, h);

    vector_compute_RHS(Nx, Ny, Cx, Cy, h, g, A);
    tdp_vector_print(100, A, stdout);
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

