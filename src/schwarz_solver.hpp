#ifndef CHP_SCHWARZ_SOLVER_H
#define CHP_SCHWARZ_SOLVER_H

#include <utility>
#include <string>

#include "cmdline.h"
#include "proc.hpp"
#include "equation.hpp"
#include "func.hpp"
#include "solver.hpp"
#include "schwarz_printer.hpp"
#include "util.h"

namespace chp {

class schwarz_solver {

public:
    enum transfer_type {
        TRANSFER_NEUMANN,
        TRANSFER_DIRICHLET,
    };

private:
    std::pair<int,int> solve_stationary(
        proc const &p, const schwarz_printer& pr, double t = 0.0);
    void solve_unstationary(proc const &p, const schwarz_printer& pr);
    bool stop_condition(int step);
    void init_mpi_type();
    void mpi_transfer_border_data_NEUMANN(proc const& p);
    void mpi_transfer_border_data_DIRICHLET(proc const& p);
    void transfer_border_data(proc const& p, transfer_type type);

protected:
    static constexpr const double EPSILON = 1e-8;

    const func& m_func;
    bool m_stationary;
    equation m_eq;
    solver_ptr m_solver;
    vec m_tmp[2];

public:
    schwarz_solver(proc& p, struct gengetopt_args_info const &opt);
    schwarz_solver(proc &p, int function, double Lx, double Ly, int Nit, double Tmax,
                   int resolutionX, int resolutionY,
                   int recouvr, const std::string& solver);

    void run(proc &p, const schwarz_printer& pr);
};

};

#endif // CHP_SCHWARZ_SOLVER_H
