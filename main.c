#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "cmdline.h"
#include "error.h"
#include "perf/perf.h"
#include "util.h"
#include "func.h"
#include "equation.h"
#include "grad.h"
#include "cblas.h"
#include "proc.h"


/* static void */
/* create_directory(char const *dirname) */
/* { */
/*     mkdir(dirname, S_IRWXU); */
/* } */

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

/* static void */
/* read_param(char const *filename, int *Nx, int *Ny, */
/*            double *Lx, double *Ly, double *D) */
/* { */
/*     FILE *f = fopen(filename, "r"); */
/*     if (f == NULL) { */
/*         perror(filename); */
/*         exit(EXIT_FAILURE); */
/*     } */
/*     fscanf(f, "%d %d %lg %lg %lg", Nx, Ny, Lx, Ly, D); */
/*     fclose(f); */
/* } */

static void
handle_opt(struct gengetopt_args_info *opt)
{
    if (opt->list_function_flag) {
        printf("Function list:\n");
        struct chp_func *list = func_list;
        int i = 0;
        while (list != NULL) {
            printf("%d: %s\n", i, list->name);
            list = list->next;
            ++ i;
        }
        exit(EXIT_SUCCESS);
    }
}

static void
chp_mpi_transfer_border_data(struct chp_proc *p, struct chp_equation *eq)
{
    int rank = p->rank, group_size = p->group_size;
    int Nx = eq->Nx, Ny = eq->Ny;
    MPI_Status st;

    if (rank < group_size-1)
        MPI_Bsend(eq->U0 + (eq->next_border_col * Nx),
                  Ny, MPI_DOUBLE, p->rank+1, rank, MPI_COMM_WORLD);
    if (rank > 0)
        MPI_Bsend(eq->U0 + (eq->prev_border_col * Nx),
                  Ny, MPI_DOUBLE, p->rank-1, rank-1, MPI_COMM_WORLD);

    if (rank < group_size-1)
        MPI_Recv(eq->right, Ny, MPI_DOUBLE, p->rank+1, rank, MPI_COMM_WORLD, &st);
    if (rank > 0)
        MPI_Recv(eq->left, Ny, MPI_DOUBLE, p->rank-1, rank-1, MPI_COMM_WORLD, &st);
}

static bool chp_stop_condition(struct chp_equation *eq, int step)
{
    int N = eq->Nx*eq->Ny;
    if (step == 0) {
        SWAP_POINTER(eq->U1, eq->U0);
        return false;
    }

    cblas_daxpy(N, -1, eq->U1, 1, eq->U0, 1);
    double n = cblas_dnrm2(N, eq->U0, 1);
    double b = cblas_dnrm2(N, eq->U1, 1);

    SWAP_POINTER(eq->U1, eq->U0);

    double v = n/b;
    double max;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (max < EPSILON) {
        SWAP_POINTER(eq->U1, eq->U0);
        return true;
    }
    return false;
}

/*stationary*/
static void
solve_equation_schwarz(struct chp_proc *p, struct gengetopt_args_info *opt)
{
    int s = opt->resolution_arg;
    int r = opt->recouvr_arg;

    struct chp_equation eq;
    struct chp_func *func;

    chp_equation_init(&eq, p->rank, p->group_size, r, s, s,
                      opt->Lx_arg, opt->Ly_arg);
    chp_equation_alloc(&eq);

    func = chp_get_func_by_id(opt->function_arg);
    chp_func_specialize_rank(func, p->rank, p->group_size);
    chp_equation_grid_init(&eq);
    chp_equation_border_init(p, &eq, func);

    if (p->rank == 1) {
        printf("Nx=%d; Ny=%d\n", eq.Nx, eq.Ny);
        printf("-1=%g\n", eq.X[0]-eq.dx);
        tdp_vector_print(eq.Nx, eq.X, stdout);
        printf("+1=%g\n", eq.X[eq.Nx-1]+eq.dx);
        printf("prev_x: %g :: next_x: %g\n",
               eq.prev_border_x, eq.next_border_x);
        printf("prev_col: %d :: next_col: %d\n",
               eq.prev_border_col, eq.next_border_col);
        printf("\n\n");
    }

    bool quit = false;
    int step = 0;
    while (!quit) {
        cblas_dcopy(eq.N, eq.rhs_f, 1, eq.rhs, 1);
        vector_compute_RHS(&eq);
        matrix_5diag_conjugate_gradient(
            eq.Nx, eq.Ny, eq.B, eq.Cx, eq.Cy, eq.rhs, eq.U1);
        chp_mpi_transfer_border_data(p, &eq);
        quit = chp_stop_condition(&eq, step);
        ++ step;
    }

    chp_output("numeric.dat", eq.Nx, eq.Ny, eq.X, eq.Y, eq.U1);
    chp_equation_free(&eq);
}

int main(int argc, char *argv[])
{
    struct gengetopt_args_info opt;

    cmdline_parser(argc, argv, &opt);
    handle_opt(&opt);

    struct chp_proc P;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &P.group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &P.rank);

    //assert( ((void)"il faut exactement 2 processus MPI", group_size == 2) );

    perf_t p1, p2;
    perf(&p1);
    solve_equation_schwarz(&P, &opt);
    perf(&p2);
    perf_diff(&p1, &p2);

    uint64_t micro, max_t;
    micro = perf_get_micro(&p2);
    MPI_Reduce(&micro, &max_t, 1, MPI_UNSIGNED_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (!P.rank)
        fprintf(stderr, "time: %lu.%06lu s\n", max_t/1000000UL, max_t%1000000UL);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
