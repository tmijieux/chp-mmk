#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

#include "error.hpp"
#include "func.hpp"
#include "solver.hpp"
#include "util.h"

using namespace chp;

std::vector<func*> func::func_list;

void zero(int const N, double const *A, double *B, double v)
{
    (void) A;
    (void) v;
    memset(B, 0, N*sizeof B[0]);
}

static void one(int const N, double const *A, double *B, double v)
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

REGISTER_FUNCTION(f=2*(x-x^2+y-y^2),
                  STATIONARY,
                  zero, zero,
                  zero, zero,
                  rhs_1, nullptr, U_1);

REGISTER_FUNCTION(f=sin(x)+cos(y),
                  STATIONARY,
                  bottom_top_2, bottom_top_2,
                  right_left_2, right_left_2,
                  rhs_2, nullptr, U_2);

REGISTER_FUNCTION(f=e^(-(x-(Lx/2)^2))*e^(-(y-(Ly/2)^2))*cos(PI/2*t),
                  UNSTATIONARY,
                  zero, zero,
                  one, one,
                  nullptr, rhs_3, nullptr);


func::func(std::string const& name, func::type type,
           rhs_t rhs, rhs_unsta_t rhs_unsta, border_t bottom,
           border_t top, border_t right, border_t left, U_t U):
    m_type(type), m_name(name), m_rhs(rhs), m_rhs_unsta(rhs_unsta),
    m_bottom(bottom), m_top(top), m_right(right), m_left(left), m_U(U)
{
    func_list.push_back(this);
}

const func& func::specialize_rank(proc const &P)
{
    if (P.rank() > 0)
        m_left = zero;

    if (P.rank() < P.size()-1)
        m_right = zero;

    return *this;
}

const func& func::get_func(unsigned idx, proc const& P)
{
    if (idx >= func_list.size())
        throw exception("Cette fonction n'existe pas");
    return func_list[idx]->specialize_rank(P);
}

void func::print_list()
{
    std::cout << "Function list:" << std::endl;
    int i = 0;
    for (auto f : func_list)
        std::cout << i++ << ": " << f->m_name << std::endl;
}
