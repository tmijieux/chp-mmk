#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/stat.h>

#include "error.hpp"
#include "util.h"
#include "cblas.h"

#include "schwarz_solver.hpp"
#include "solver.hpp"

using namespace chp;

static void
create_directory(char const *dirname)
{
    mkdir(dirname, S_IRWXU);
}

static void
output(char const *filename, int const Nx, int const Ny,
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

void schwarz_solver::init_mpi_type()
{
    MPI_Type_vector(_eq.Ny, 1, _eq.Nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
}

void schwarz_solver::mpi_transfer_border_data_NEUMANN(proc const& p)
{
    int rank = p.rank(), group_size = p.size();
    int Ny = _eq.Ny, Nx = _eq.Nx;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    double *tmp_right = &_tmp[0][0];
    double *tmp_left = &_tmp[1][0];

    cblas_dcopy(Ny, _eq.U1+_eq.next_border_col2, Nx, tmp_right, 1);
    cblas_daxpy(Ny, -1.0, _eq.U1+_eq.next_border_col, Nx, tmp_right, 1);
    cblas_dcopy(Ny, _eq.U1+_eq.prev_border_col2, Nx, tmp_left, 1);
    cblas_daxpy(Ny, -1.0, _eq.U1+_eq.prev_border_col, Nx, tmp_left, 1);

    if (rank < group_size-1)
        MPI_Isend(tmp_right, Ny, MPI_DOUBLE, p.rank()+1,
                  rank, MPI_COMM_WORLD, &r[0]);

    if (rank > 0)
        MPI_Isend(tmp_left, Ny, MPI_DOUBLE, p.rank()-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(_eq.right, Ny, MPI_DOUBLE,
                  p.rank()+1, group_size+rank, MPI_COMM_WORLD, &r[2]);

    if (rank > 0)
        MPI_Irecv(_eq.left, Ny, MPI_DOUBLE,
                  p.rank()-1, rank-1, MPI_COMM_WORLD, &r[3]);

    MPI_Waitall(4, r, st);

    cblas_daxpy(Ny, 1.0, _eq.U1+Ny-1, Nx, _eq.right, 1);
    cblas_daxpy(Ny, 1.0, _eq.U1, Nx, _eq.left, 1);
}

void schwarz_solver::mpi_transfer_border_data_DIRICHLET(proc const& p)
{
    int rank = p.rank(), group_size = p.size();
    int Ny = _eq.Ny;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    if (rank < group_size-1)
        MPI_Isend(_eq.U1 + _eq.next_border_col,
                  1, column_type, p.rank()+1, rank, MPI_COMM_WORLD, &r[0]);
    if (rank > 0)
        MPI_Isend(_eq.U1 + _eq.prev_border_col,
                  1, column_type, p.rank()-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(_eq.right, Ny, MPI_DOUBLE,
                  p.rank()+1, group_size+rank, MPI_COMM_WORLD, &r[2]);
    if (rank > 0)
        MPI_Irecv(_eq.left, Ny, MPI_DOUBLE,
                  p.rank()-1, rank-1, MPI_COMM_WORLD, &r[3]);
    MPI_Waitall(4, r, st);
}

void schwarz_solver::transfer_border_data(proc const& p, transfer_type type)
{
    if (type == TRANSFER_NEUMANN)
        mpi_transfer_border_data_NEUMANN(p);
    else if  (type == TRANSFER_DIRICHLET)
        mpi_transfer_border_data_DIRICHLET(p);
}

bool schwarz_solver::stop_condition(int step)
{
    int N = _eq.Nx*_eq.Ny;
    if (step == 0) {
        std::swap(_eq.U1, _eq.U0);
        return false;
    }

    // U0 = ancienne solution
    //   (calculé au tour précédant ou vecteur nul pour le tour 0)
    // U1 = solution fraichement calculé
    cblas_daxpy(N, -1, _eq.U1, 1, _eq.U0, 1); // U0 = U0 - U1
    double n = cblas_dnrm2(N, _eq.U0, 1);
    double b = cblas_dnrm2(N, _eq.U1, 1);

    std::swap(_eq.U1, _eq.U0);
    // U0 pointe desormais vers la solution "fraichement calculé"
    // et U1 vers la différence des anciens U0-U1

    double v = n/b;
    double max;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (max < EPSILON) {
        std::swap(_eq.U1, _eq.U0);
        return true;
    }
    return false;
}

void schwarz_solver::solve_stationary(proc const &p)
{
    bool quit = false;
    int schwarz_step = 0;
    int total_step = 0;

    _eq.rhs_init(_func, 0.0);

    while (!quit) {
        _eq.apply_border_cond_RHS();
        int s = _S->solve(_eq.rhs, _eq.U0, _eq.U1);
        if (!p.rank())
            printf("schwarz_step(%d): solver_step = %d\n", schwarz_step, s);

        total_step += s;
        transfer_border_data(p, TRANSFER_DIRICHLET);
        quit = stop_condition(schwarz_step);
        ++ schwarz_step;
    }
    if (!p.rank())
        printf("# schwarz_step = %d\n# total_step = %d\n", schwarz_step, total_step);

    char filename[100];
    snprintf(filename, 100, "numeric.dat.%d", p.rank());
    output(filename, _eq.Nx, _eq.Ny, _eq.X, _eq.Y, _eq.U1);
}

void schwarz_solver::solve_unstationary(proc const& p)
{
    const int print_step = 100;
    double t = 0.0;
    int total_step = 0;
    int total_schwarz_step = 0;

    create_directory("sol");

    for (int i = 0; i < _eq.Nit; ++i) {
        t += _eq.dt;
        _eq.rhs_init(_func, t);

        int schwarz_step = 0;
        bool quit = false;
        while (!quit) {
            _eq.apply_border_cond_RHS();
            total_step += _S->solve(_eq.rhs, _eq.U0, _eq.U1);
            transfer_border_data(p, TRANSFER_DIRICHLET);
            quit = stop_condition(schwarz_step);
            ++ schwarz_step;
        }
        total_schwarz_step += schwarz_step;

        if ((i % print_step) == 0) {
            if (!p.rank())
                printf("# time_step(%d): schwarz_step = %d\n", i, schwarz_step);

            char filename[100];
            snprintf(filename, 100,
                     "sol/sol%d.dat.%d", i/print_step + 1, p.rank());
            output(filename, _eq.Nx, _eq.Ny, _eq.X, _eq.Y, _eq.U1);
        }
    }

    if (!p.rank())
        printf("# total_schwarz_step = %d\n# total_step = %d\n",
               total_schwarz_step, total_step);
}

schwarz_solver::schwarz_solver(proc &p, struct gengetopt_args_info *opt):
    _func(func::get_func(opt->function_arg, p)),
    _stationary(_func._type == func::STATIONARY),
    _eq(p, opt, _stationary),
    _S(solver::make_solver(_eq, opt->solver_arg))
{
    if (_func._type != func::STATIONARY && _func._type != func::UNSTATIONARY)
        throw exception("invalid method type");

    _eq.border_init(_func);

    init_mpi_type();

    for (int i = 0; i < 2; ++i)
        _tmp[i].resize(_eq.Ny);
}

void schwarz_solver::run(proc &p)
{
    if (_stationary)
        solve_stationary(p);
    else
        solve_unstationary(p);
}
