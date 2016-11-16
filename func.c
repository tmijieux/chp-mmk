#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "func.h"

unsigned method_list_length = 0;
struct projchp_method *method_list = NULL;

static void
zero(int const N, double const *A, double *B, double param)
{
    (void) A;
    (void) param;
    memset(B, 0, 2*N*sizeof B[0]);
}

static void
one(int const N, double const *A, double *B, double param)
{
    (void) A;
    (void) param;
    for (int i = 0; i < 2*N; ++i)
        B[i] = 1;
}

static void
f1(int const Nx, int const Ny,
   double const *X, double const *Y,
   double *F)
{
    for (int j = 0; j < Ny; ++j) {
        double yy = Y[j] - SQUARE(Y[j]);
        for (int i = 0; i < Nx; ++i)
            F[Nx*j+i] = 2*(yy + X[i] - SQUARE(X[i]));
    }
}

static void
mU1(int const Nx, int const Ny,
    double const *X, double const *Y,
    double const Lx, double const Ly, double *U)
{
    for (int j = 0; j < Ny; ++j) {
        double yy = Y[j] * (Ly-Y[j]);
        for (int i = 0; i < Nx; ++i)
            U[Nx*j+i] = X[i] * (Lx-X[i]) * yy;
    }
}

static void
g2(int const Nx, double const *X, double *G, double Lx_min, double Ly)
{
    double c = cos(Lx_min);
    for (int i = 0; i < Nx; ++i)
        G[i] = sin(X[i]) + c;

    c = cos(Ly);
    for (int i = 0; i < Nx; ++i)
        G[Nx+i] = sin(X[i]) + c;
}

static void
h2(int const Ny, double const *Y, double *G, double Lx_min, double Lx_max)
{
    double s = sin(Lx_min);
    for (int i = 0; i < Ny; ++i)
        G[i] = cos(Y[i]) + s;

    double c = sin(Lx_max);
    for (int i = 0; i < Ny; ++i)
        G[Ny+i] = cos(Y[i]) + c;
}

static void f2(int const Nx, int const Ny,
               double const *X, double const *Y,
               double *F)
{
    for (int j = 0; j < Ny; ++j) {
        double s = cos(Y[j]);
        for (int i = 0; i < Nx; ++i)
            F[Nx*j+i] = sin(X[i]) + s;
    }
}

static void
mU2(int const Nx, int const Ny,
    double const *X, double const *Y,
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
f3(int const Nx, int const Ny,
   double const *X, double const *Y,
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
                  zero, one, NULL, f3, NULL);

REGISTER_FUNCTION(f=sin(x)+cos(y),
                  STATIONARY,
                  g2, h2, f2, NULL, mU2);

REGISTER_FUNCTION(f=2*(x-x^2+y-y^2),
                  STATIONARY,
                  zero, zero, f1, NULL, mU1);
