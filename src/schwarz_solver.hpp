#ifndef CHP_SCHWARZ_SOLVER_H
#define CHP_SCHWARZ_SOLVER_H

#include "cmdline.h"

#include "proc.hpp"
#include "equation.hpp"
#include "func.hpp"
#include "solver.hpp"

namespace chp {

class schwarz_solver {

public:
    enum transfer_type {
        TRANSFER_NEUMANN,
        TRANSFER_DIRICHLET,
    };

private:
    void solve_stationary(proc const &p);
    void solve_unstationary(proc const &p);
    bool stop_condition(int step);
    void init_mpi_type();
    void mpi_transfer_border_data_NEUMANN(proc const& p);
    void mpi_transfer_border_data_DIRICHLET(proc const& p);
    void transfer_border_data(proc const& p, transfer_type type);


protected:
    static constexpr const double EPSILON = 1e-8;

    const func& _func;
    bool _stationary;
    equation _eq;
    solver *_S;


    std::vector<double> _tmp[2];

public:
    schwarz_solver(proc& p, struct gengetopt_args_info *opt);
    void run(proc &p);
};

};

#endif // CHP_SCHWARZ_SOLVER_H
