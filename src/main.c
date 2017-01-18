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

typedef enum chp_transfer_type_ chp_transfer_type;
enum chp_transfer_type_ {
    CHP_TRANSFER_NEUMANN,
    CHP_TRANSFER_DIRICHLET,
};

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

MPI_Datatype column_type;

static void
chp_mpi_init_type(struct chp_equation *eq)
{
    MPI_Type_vector(eq->Ny, 1, eq->Nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
}

static void
chp_mpi_transfer_border_data_NEUMANN(
    struct chp_proc *p, struct chp_equation *eq)
{
    int rank = p->rank, group_size = p->group_size;
    int Ny = eq->Ny, Nx = eq->Nx;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

    double *tmp_right, *tmp_left;
    MALLOC_ARRAY(tmp_right, eq->Ny);
    MALLOC_ARRAY(tmp_left, eq->Ny);

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

    free(tmp_right);
    free(tmp_left);
}

static void
chp_mpi_transfer_border_data_DIRICHLET(
    struct chp_proc *p, struct chp_equation *eq)
{
    int rank = p->rank, group_size = p->group_size;
    int Ny = eq->Ny;
    MPI_Request r[4] = {
        MPI_REQUEST_NULL, MPI_REQUEST_NULL,
        MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status st[4];
    memset(st, 0, sizeof st);

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
chp_mpi_transfer_border_data(
    struct chp_proc *p, struct chp_equation *eq,
    chp_transfer_type ttype)
{
    if (ttype == CHP_TRANSFER_NEUMANN)
        chp_mpi_transfer_border_data_NEUMANN(p, eq);
    else if  (ttype == CHP_TRANSFER_DIRICHLET)
        chp_mpi_transfer_border_data_DIRICHLET(p, eq);
}

static bool
chp_stop_condition(struct chp_equation *eq, struct chp_proc *p, int step)
{
    int N = eq->Nx*eq->Ny;
    if (step == 0) {
        SWAP_POINTER(eq->U1, eq->U0);
        return false;
    }

    // U0 = ancienne solution (calculé au tour précédant
    //                 ou vecteur nul pour le tour 0)
    // U1 = solution fraichement calculé
    cblas_daxpy(N, -1, eq->U1, 1, eq->U0, 1); // U0 = U0 - U1
    double n = cblas_dnrm2(N, eq->U0, 1);
    double b = cblas_dnrm2(N, eq->U1, 1);
    // cracra

    SWAP_POINTER(eq->U1, eq->U0);
    // U0 pointe desormais vers la solution "fraichement calculé"
    // et U1 vers la différence des anciens U0-U1

    double v = n/b;
    double max;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i< p->group_size; ++i){
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == p->rank)
            printf("stop value(%d): ||(U0-U1)|| / ||U0|| = %g\n", p->rank, v);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (p->rank == 0)
        printf("stop value: ||(U0-U1)|| / ||U0|| = %g\n\n", max);
    MPI_Barrier(MPI_COMM_WORLD);
    if (p->rank == 0)
        puts("");
    MPI_Barrier(MPI_COMM_WORLD);

    if (max < EPSILON) {
        SWAP_POINTER(eq->U1, eq->U0);
        return true;
    }
    return false;
}

void print_debug_info_1(struct chp_proc *p, struct chp_equation *eq)
{
    if (p->rank == 1) {
        printf("Nx=%d; Ny=%d\n", eq->Nx, eq->Ny);
        printf("-1=%g\n", eq->X[0]-eq->dx);
        tdp_vector_print(eq->Nx, eq->X, stdout);
        printf("+1=%g\n", eq->X[eq->Nx-1]+eq->dx);
        printf("prev_x: %g :: next_x: %g\n",
               eq->prev_border_x, eq->next_border_x);
        printf("prev_col: %d :: next_col: %d\n",
               eq->prev_border_col, eq->next_border_col);
        printf("\n\n");
    }
}

void print_debug_info_2(
    struct chp_proc *p, struct chp_equation *eq, int step)
{
    if (p->rank == 1) {
        printf("left step %d:\n", step);
        tdp_vector_print(eq->Ny, eq->left, stdout);

        printf("right step %d:\n", step);
        tdp_vector_print(eq->Ny, eq->right, stdout);
    }
}

/*stationary*/
static void
solve_equation_schwarz_stationary(
    struct chp_proc *p, struct chp_equation *eq, struct chp_func *func, char* solver)
{
    chp_equation_border_init(eq, func);
    chp_equation_rhs_init(eq, func, 0.0);

    bool quit = false;
    int step = 0;
    while (!quit) {
        cblas_dcopy(eq->N, eq->rhs_f, 1, eq->rhs, 1);
        vector_compute_RHS(eq);

        if(!strcmp(solver,"CG"))
            matrix_5diag_conjugate_gradient(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
	else if(!strcmp(solver,"GS"))
            matrix_5diag_gauss_seidel(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
	else if(!strcmp(solver,"J"))
            matrix_5diag_jacobi(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
	else{
            fprintf(stderr, "Solveur invalide: '%s'\n", solver);
            exit(EXIT_FAILURE);
	}

        //print_debug_info_2(p, &eq, step);

        chp_mpi_transfer_border_data(p, eq, CHP_TRANSFER_DIRICHLET);
        quit = chp_stop_condition(eq, p, step);
        ++ step;
    }
    if (!p->rank)
        printf("#step=%d\n", step);

    char filename[100];
    snprintf(filename, 100, "numeric.dat.%d", p->rank);
    chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
}

/*unstationary*/
static void
solve_equation_schwarz_unstationary(
    struct chp_proc *p, struct chp_equation *eq, struct chp_func *func, char* solver)
{
    create_directory("sol");
    chp_equation_border_init(eq, func);
    double const Tmax = 10.0;
    int const Nit = 2000;
    double dt = Tmax / Nit;
    const int print_step = 1;
    eq->B += 1.0 / dt;
    eq->dt = dt;

    double t = 0.0;
    for (int i = 0; i < Nit; ++i) {
        t += dt;
        chp_equation_rhs_init(eq, func, t);
        int step = 0;
        bool quit = false;
        while (!quit) {
            cblas_dcopy(eq->N, eq->rhs_f, 1, eq->rhs, 1);
            vector_compute_RHS(eq);

            if(!strcmp(solver,"CG"))
                matrix_5diag_conjugate_gradient(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
            else if(!strcmp(solver,"GS"))
                matrix_5diag_gauss_seidel(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
            else if(!strcmp(solver,"J"))
                matrix_5diag_jacobi(eq->Nx, eq->Ny, eq->B, eq->Cx, eq->Cy, eq->rhs, eq->U0, eq->U1);
            else{
                fprintf(stderr, "Solveur invalide: '%s'\n", solver);
                exit(EXIT_FAILURE);
            }

            chp_mpi_transfer_border_data(p, eq, CHP_TRANSFER_DIRICHLET);
            quit = chp_stop_condition(eq, p, step);
            ++ step;
        }

        if ((i % print_step) == 0) {
            if (!p->rank)
                printf("#step = %d\n", step);

            char filename[100];
            snprintf(filename, 100,
                     "sol/sol%d.dat.%d", i/print_step + 1, p->rank);
            chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
            /* printf("i=%d\n", i); */
        }
    }
}

static void solve_equation_schwarz(
    struct chp_proc *p, struct gengetopt_args_info *opt)
{
    int NX = opt->resolutionX_arg;
    int NY = opt->resolutionY_arg;
    int r = opt->recouvr_arg;
    char* solver = opt->solver_arg;

    struct chp_equation eq;
    struct chp_func *func;

    chp_equation_init(&eq, p->rank, p->group_size, r, NX, NY,
                      opt->Lx_arg, opt->Ly_arg);
    chp_mpi_init_type(&eq);
    chp_equation_alloc(&eq);

    func = chp_get_func_by_id(opt->function_arg);
    chp_func_specialize_rank(func, p->rank, p->group_size);
    chp_equation_grid_init(&eq);

    if (func->type == CHP_UNSTATIONARY)
        solve_equation_schwarz_unstationary(p, &eq, func, solver);
    else if (func->type == CHP_STATIONARY)
        solve_equation_schwarz_stationary(p, &eq, func, solver);
    else
        projchp_error("invalid method type");

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

    //assert( ((void)"il faut exactement 2 processus MPI",
    //         group_size == 2) );

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
        fprintf(stderr, "# %dx%d: %lu.%06lu s\n",
                opt.resolutionX_arg, opt.resolutionY_arg,
                max_t/1000000UL, max_t%1000000UL);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
