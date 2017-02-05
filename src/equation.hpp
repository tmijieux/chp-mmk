#ifndef CHP_EQUATION_H
#define CHP_EQUATION_H

#include "func.hpp"
#include "proc.hpp"
#include "cmdline.h"
#include "util.h"

namespace chp {

class equation {

private:

    void init_grid();
    void init_prop(proc const& P, int, int, int, bool);
    void alloc();

public:
    int Nx, Ny, N;
    double Lx_min, Lx_max;
    double Ly_min, Ly_max;
    double Lx, Ly;

    int Nit;
    double Tmax;
    double dx, dy, dt;
    double D, B, Cx, Cy;

    vec X, Y;
    vec top, bottom, left, right;
    vec rhs, rhs_f;
    vec U0, U1;

    double next_border_x, prev_border_x;
    int next_border_col, prev_border_col;
    int next_border_col2, prev_border_col2;

    equation(proc const& P, struct gengetopt_args_info const& opt, bool stationary);
    void border_init(func const &f);
    void rhs_init(func const &f, double t);
    void apply_border_cond_RHS();
};

};


#endif // CHP_EQUATION_H
