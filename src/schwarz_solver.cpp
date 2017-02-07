#include <cstdlib>
#include <cstring>
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

MPI_Datatype column_type;

void schwarz_solver::init_mpi_type()
{
    MPI_Type_vector(m_eq.Ny, 1, m_eq.Nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
}

void schwarz_solver::mpi_transfer_border_data_NEUMANN(const proc& p)
{
    int rank = p.rank(), group_size = p.size();
    int Ny = m_eq.Ny, Nx = m_eq.Nx;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    double *tmp_right = m_tmp[0];
    double *tmp_left = m_tmp[1];

    cblas_dcopy(Ny, m_eq.U1+m_eq.next_border_col2, Nx, tmp_right, 1);
    cblas_daxpy(Ny, -1.0, m_eq.U1+m_eq.next_border_col, Nx, tmp_right, 1);
    cblas_dcopy(Ny, m_eq.U1+m_eq.prev_border_col2, Nx, tmp_left, 1);
    cblas_daxpy(Ny, -1.0, m_eq.U1+m_eq.prev_border_col, Nx, tmp_left, 1);

    if (rank < group_size-1)
        MPI_Isend(tmp_right, Ny, MPI_DOUBLE, p.rank()+1,
                  rank, MPI_COMM_WORLD, &r[0]);

    if (rank > 0)
        MPI_Isend(tmp_left, Ny, MPI_DOUBLE, p.rank()-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(m_eq.right, Ny, MPI_DOUBLE,
                  p.rank()+1, group_size+rank, MPI_COMM_WORLD, &r[2]);

    if (rank > 0)
        MPI_Irecv(m_eq.left, Ny, MPI_DOUBLE,
                  p.rank()-1, rank-1, MPI_COMM_WORLD, &r[3]);

    MPI_Waitall(4, r, st);

    cblas_daxpy(Ny, 1.0, m_eq.U1+Ny-1, Nx, m_eq.right, 1);
    cblas_daxpy(Ny, 1.0, m_eq.U1, Nx, m_eq.left, 1);
}

void schwarz_solver::mpi_transfer_border_data_DIRICHLET(const proc& p)
{
    int rank = p.rank(), group_size = p.size();
    int Ny = m_eq.Ny;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    if (rank < group_size-1)
        MPI_Isend(m_eq.U1 + m_eq.next_border_col,
                  1, column_type, p.rank()+1, rank, MPI_COMM_WORLD, &r[0]);
    if (rank > 0)
        MPI_Isend(m_eq.U1 + m_eq.prev_border_col,
                  1, column_type, p.rank()-1,
                  group_size+rank-1, MPI_COMM_WORLD, &r[1]);

    if (rank < group_size-1)
        MPI_Irecv(m_eq.right, Ny, MPI_DOUBLE,
                  p.rank()+1, group_size+rank, MPI_COMM_WORLD, &r[2]);
    if (rank > 0)
        MPI_Irecv(m_eq.left, Ny, MPI_DOUBLE,
                  p.rank()-1, rank-1, MPI_COMM_WORLD, &r[3]);
    MPI_Waitall(4, r, st);
}

void schwarz_solver::transfer_border_data(const proc& p, transfer_type type)
{
    if (type == TRANSFER_NEUMANN)
        mpi_transfer_border_data_NEUMANN(p);
    else if  (type == TRANSFER_DIRICHLET)
        mpi_transfer_border_data_DIRICHLET(p);
}

bool schwarz_solver::stop_condition(int step)
{
    int N = m_eq.N;
    if (step == 0) {
        std::swap(m_eq.U1, m_eq.U0);
        return false;
    }

    // U0 = ancienne solution
    //   (calculé au tour précédant ou vecteur nul pour le tour 0)
    // U1 = solution fraichement calculé
    cblas_daxpy(N, -1, m_eq.U1, 1, m_eq.U0, 1); // U0 = U0 - U1
    double n = cblas_dnrm2(N, m_eq.U0, 1);
    double b = cblas_dnrm2(N, m_eq.U1, 1);

    std::swap(m_eq.U1, m_eq.U0);
    // U0 pointe desormais vers la solution "fraichement calculé"
    // et U1 vers la différence des anciens U0-U1

    double v = n/b;
    double max;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (max < EPSILON) {
        std::swap(m_eq.U1, m_eq.U0);
        return true;
    }
    return false;
}

std::pair<int,int>
schwarz_solver::solve_stationary(const proc &p, const schwarz_printer& pr, double t)
{
    bool quit = false;
    int schwarz_step = 0;
    int total_step = 0;

    m_eq.rhs_init(m_func, t);

    while (!quit) {
        m_eq.apply_border_cond_RHS();
        int step = m_solver->solve(m_eq.rhs, m_eq.U0, m_eq.U1);
        pr.print_sta_step(schwarz_step, step);
        total_step += step;
        transfer_border_data(p, TRANSFER_DIRICHLET);
        quit = stop_condition(schwarz_step);
        ++ schwarz_step;
    }
    pr.print_sta_fin(schwarz_step, total_step, m_eq);
    return std::make_pair(schwarz_step, total_step);
}

void schwarz_solver::solve_unstationary(const proc& p, const schwarz_printer& pr)
{
    double t = 0.0;
    int step = 0, sstep = 0;
    schwarz_printer silent_printer(p);

    create_directory("sol");

    for (int i = 0; i < m_eq.Nit; ++i) {
        t += m_eq.dt;
        int s, ss;
        std::tie(s, ss) = solve_stationary(p, silent_printer, t);
        step += s; sstep += ss;
        pr.print_unsta_step(i, ss, m_eq);
    }
    pr.print_unsta_fin(step, sstep);
}


schwarz_solver::
schwarz_solver( proc &p, int function,
                int resolutionX, int resolutionY,
                int recouvr, const string& solver  ):
    
    m_func(func::get_func(function, p)),
    m_stationary(m_func.m_type == func::STATIONARY),
    m_eq(p, resolutionX, resolutionY, recouvr, m_stationary),
    m_solver(solver::make_solver(m_eq, solver))
    
{
    if (m_func.m_type != func::STATIONARY && m_func.m_type != func::UNSTATIONARY)
        throw exception("invalid method type");

    m_eq.border_init(m_func);
    init_mpi_type();

    for (int i = 0; i < 2; ++i)
        m_tmp[i].resize(m_eq.Ny);
}


schwarz_solver::schwarz_solver(proc &p, struct gengetopt_args_info const &opt):
    schwarz_solver(
        p, opt.function_arg,
        opt.resolution_X_arg, opt.resolutionY_arg,
        opt.recouvr_arg, opt.solver_arg
    )
{
}

void schwarz_solver::run(proc &p, const schwarz_printer& pr)
{
    if (m_stationary)
        solve_stationary(p, pr);
    else
        solve_unstationary(p, pr);
}
