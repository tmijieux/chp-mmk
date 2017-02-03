#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/stat.h>

#include "util.h"
#include "cblas.h"

#include "schwarz_solver.h"
#include "solver.h"

static void
create_directory(char const *dirname)
{
    mkdir(dirname, S_IRWXU);
}

static void
chp_output(char const *filename, int const Nx, int const Ny,
           double const *X, double const *Y, double const *U)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j)
            fprintf(f, "%g %g %g\n", X[i], Y[j], U[Nx*j+i]);
        fprintf(f, "\n");
    }
    fclose(f);
}

MPI_Datatype column_type;

static void chp_mpi_init_type(const chp_equation *eq)
{
    MPI_Type_vector(eq->Ny, 1, eq->Nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
}

static void
chp_mpi_transfer_border_data_NEUMANN(chp_proc *p, chp_equation *eq, double *restrict tmp[])
{
    int rank = p->rank, group_size = p->group_size;
    int Ny = eq->Ny, Nx = eq->Nx;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    double *tmp_right = tmp[0];
    double *tmp_left = tmp[1];

    cblas_dcopy(Ny, eq->U1+eq->next_border_col2, Nx, tmp_right, 1);
    cblas_daxpy(Ny, -1.0, eq->U1+eq->next_border_col, Nx, tmp_right, 1);
    cblas_dcopy(Ny, eq->U1+eq->prev_border_col2, Nx, tmp_left, 1);
    cblas_daxpy(Ny, -1.0, eq->U1+eq->prev_border_col, Nx, tmp_left, 1);

    if (rank < group_size-1)
        MPI_Isend(tmp_right, Ny, MPI_DOUBLE, p->rank+1,
                  rank, MPI_COMM_WORLD, &r[0]);

    if (rank > 0)
        MPI_Isend(tmp_left, Ny, MPI_DOUBLE, p->rank-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(eq->right, Ny, MPI_DOUBLE,
                  p->rank+1, group_size+rank, MPI_COMM_WORLD, &r[2]);

    if (rank > 0)
        MPI_Irecv(eq->left, Ny, MPI_DOUBLE,
                  p->rank-1, rank-1, MPI_COMM_WORLD, &r[3]);

    MPI_Waitall(4, r, st);

    cblas_daxpy(Ny, 1.0, eq->U1+Ny-1, Nx, eq->right, 1);
    cblas_daxpy(Ny, 1.0, eq->U1, Nx, eq->left, 1);
}

static void
chp_mpi_transfer_border_data_DIRICHLET(chp_proc *p, chp_equation *eq, double * restrict tmp[])
{
    int rank = p->rank, group_size = p->group_size;
    int Ny = eq->Ny;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);
    (void) tmp;

    if (rank < group_size-1)
        MPI_Isend(eq->U1 + eq->next_border_col,
                  1, column_type, p->rank+1, rank, MPI_COMM_WORLD, &r[0]);
    if (rank > 0)
        MPI_Isend(eq->U1 + eq->prev_border_col,
                  1, column_type, p->rank-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(eq->right, Ny, MPI_DOUBLE,
                  p->rank+1, group_size+rank, MPI_COMM_WORLD, &r[2]);
    if (rank > 0)
        MPI_Irecv(eq->left, Ny, MPI_DOUBLE,
                  p->rank-1, rank-1, MPI_COMM_WORLD, &r[3]);
    MPI_Waitall(4, r, st);
}

static void
chp_mpi_transfer_border_data(chp_proc *p, chp_equation *eq,
                             chp_transfer_type type, double *restrict tmp[])
{
    if (type == CHP_TRANSFER_NEUMANN)
        chp_mpi_transfer_border_data_NEUMANN(p, eq, tmp);
    else if  (type == CHP_TRANSFER_DIRICHLET)
        chp_mpi_transfer_border_data_DIRICHLET(p, eq, tmp);
}

static bool
chp_stop_condition(chp_equation *eq, int step)
{
    int N = eq->Nx*eq->Ny;
    if (step == 0) {
        SWAP_POINTER(eq->U1, eq->U0);
        return false;
    }

    // U0 = ancienne solution
    //   (calculé au tour précédant ou vecteur nul pour le tour 0)
    // U1 = solution fraichement calculé
    cblas_daxpy(N, -1, eq->U1, 1, eq->U0, 1); // U0 = U0 - U1
    double n = cblas_dnrm2(N, eq->U0, 1);
    double b = cblas_dnrm2(N, eq->U1, 1);

    SWAP_POINTER(eq->U1, eq->U0);
    // U0 pointe desormais vers la solution "fraichement calculé"
    // et U1 vers la différence des anciens U0-U1

    double v = n/b;
    double max;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (max < SCHWARZ_EPSILON) {
        SWAP_POINTER(eq->U1, eq->U0);
        return true;
    }
    return false;
}

static void
solve_stationary(chp_proc *p, chp_schwarz_solver *sS, double t, bool output)
{
    chp_equation *eq = &sS->eq;
    chp_func *func = &sS->func;
    chp_solver *S = &sS->S;

    chp_equation_rhs_init(eq, func, t);
    
    bool quit = false;    
    int schwarz_step = 0;
    while (!quit) {
        chp_equation_apply_border_cond_RHS(eq);
        int step = chp_solver_run(S, eq->rhs, eq->U0, eq->U1);
        chp_mpi_transfer_border_data(p, eq, CHP_TRANSFER_DIRICHLET, sS->tmp);
        quit = chp_stop_condition(eq, schwarz_step);
        ++ schwarz_step;
    }

    if (output) {
        char filename[100];
        snprintf(filename, 100, "numeric.dat.%d", p->rank);
        chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
    }
}

static void
solve_unstationary(chp_proc *p, chp_schwarz_solver *sS, bool output)
{
    const int print_step = 100;
    double t = 0.0;
    chp_equation *eq = &sS->eq;
    create_directory("sol");

    for (int i = 0; i < eq->Nit; ++i) {
        t += eq->dt;
        solve_stationary(p, sS, t, false);

        if (output) {
            if ((i % print_step) == 0) {
                char filename[100];
                snprintf(filename, 100,
                         "sol/sol%d.dat.%d", i/print_step + 1, p->rank);
                chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
            }
        }
    }
    if (output) {
        
    }
}

void chp_schwarz_solver_init(
    chp_schwarz_solver *S, chp_proc *p, struct gengetopt_args_info *opt)
{
    memset(S, 0, sizeof*S);

    chp_func_init(&S->func, opt->function_arg, p);
    S->stationary = (S->func.type == CHP_STATIONARY);

    if (S->func.type != CHP_STATIONARY && S->func.type != CHP_UNSTATIONARY)
        chp_error("invalid method type");

    chp_equation_init(&S->eq, p, opt, S->stationary);
    chp_equation_border_init(&S->eq, &S->func);
    chp_mpi_init_type(&S->eq);
    chp_solver_init(&S->S, &S->eq, opt->solver_arg);

    for (int i = 0; i < 2; ++i)
        S->tmp[i] = tdp_vector_new(S->eq.Ny);
}

void chp_schwarz_solver_run(chp_schwarz_solver *S, chp_proc *p, bool voutput)
{
    if (S->stationary)
        solve_stationary(p, S, 0.0, voutput);
    else
        solve_unstationary(p, S, voutput);
}

void chp_schwarz_solver_free(chp_schwarz_solver *S)
{
    chp_equation_free(&S->eq);
    chp_solver_free(&S->S);
    for (int i = 0; i < 2; ++i)
        free(S->tmp[i]);

    memset(S, 0, sizeof*S);
}
