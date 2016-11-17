#include <stdio.h>
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

static void
read_param(char const *filename, int *Nx, int *Ny,
           double *Lx, double *Ly, double *D)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        perror(filename);
        exit(EXIT_FAILURE);
    }
    fscanf(f, "%d %d %lg %lg %lg", Nx, Ny, Lx, Ly, D);
    fclose(f);
}

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
    
    /* if (p->rank < p->group_size - 1) { */
    /*     MPI_Send(eq->U, Ny, MPI_DOUBLE, p->rank+1, 0, MPI_COMM_WORLD); */
    /*     MPI_Recv(eq->U, Ny, MPI_DOUBLE, p->rank+1, 0, MPI_COMM_WORLD); */
    /* } */
    /* if (p->rank > 0) { */
    /*     MPI_Recv(eq->U, Ny, MPI_DOUBLE, p->rank-1, 0, MPI_COMM_WORLD); */
    /*     MPI_Send(eq->U, Ny, MPI_DOUBLE, p->rank-1, 0, MPI_COMM_WORLD); */
    /* } */
}

/*stationary*/
static void
solve_equation_schwarz(struct chp_proc *p, struct gengetopt_args_info *opt)
{
    int s = opt->resolution_arg;
    int r = opt->recouvr_arg;

    struct chp_equation eq;
    struct chp_func *func;

    chp_equation_init(&eq, p->rank, p->group_size, r, s);
    chp_equation_alloc(&eq);
    func = chp_get_func_by_id(opt->function_arg);
    chp_func_specialize_rank(func, p->rank, p->group_size);

    chp_equation_grid_init(&eq);
    chp_equation_border_init(p, &eq, func);

    if (p->rank == 1) {
        printf("Nx=%d; Ny=%d\n", eq.Nx, eq.Ny);
        tdp_vector_print(eq.Nx, eq.X, stdout);
        printf("\n\n");
        //tdp_vector_print(eq.Ny, eq.Y, stdout);
    }        

    #define NOT_CONVERGED 0
    while ( NOT_CONVERGED) {
        matrix_5diag_conjugate_gradient(
            eq.Nx, eq.Ny, eq.B, eq.Cx, eq.Cy, eq.rhs, eq.U0);
        chp_mpi_transfer_border_data(p, &eq);
    }
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
