#ifndef CHP_SCHWARZ_SOLVER_H
#define CHP_SCHWARZ_SOLVER_H

#include <stdbool.h>

typedef enum chp_transfer_type_ chp_transfer_type;
typedef struct chp_schwarz_solver_ chp_schwarz_solver;

#include "proc.h"
#include "equation.h"
#include "func.h"
#include "cmdline.h"
#include "solver.h"


#define SCHWARZ_EPSILON   1e-8

enum chp_transfer_type_ {
    CHP_TRANSFER_NEUMANN,
    CHP_TRANSFER_DIRICHLET,
};

struct chp_schwarz_solver_ {
    bool stationary;
    chp_solver S;
    chp_equation eq;
    chp_func func;

};

void chp_schwarz_solver_init(
    chp_schwarz_solver *S, chp_proc *p, struct gengetopt_args_info *opt);
void chp_schwarz_solver_run(chp_schwarz_solver *S, chp_proc *p);
void chp_schwarz_solver_free(chp_schwarz_solver *S);

#endif // CHP_SCHWARZ_SOLVER_H
