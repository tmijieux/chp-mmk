#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <cblas.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "cmdline.h"
#include "error.h"
#include "perf/perf.h"
#include "util.h"
#include "func.h"
#include "grad.h"

static void
create_directory(char const *dirname)
{
    mkdir(dirname, S_IRWXU);
}

static void
projchp_method_call(struct projchp_method *method,
                    int const Nx, int const Ny,
                    double *Uexact, double *X, double *Y,
                    double *f, double *g, double *h,
                    double const Lx, double const Ly)
{
    method->g(Nx, X, g, Ly);
    method->h(Ny, Y, h, Lx);
    method->f(Nx, Ny, X, Y, f);
    method->mU(Nx, Ny, X, Y, Lx, Ly, Uexact);
}

static void
projchp_grid_init(double *X, double *Y,
                  double dx, double dy,
                  int const Nx, int const Ny)
{
    for (int i = 0; i < Nx; ++i)
        X[i] = i*dx;
    for (int i = 0; i < Ny; ++i)
        Y[i] = i*dy;
}

static void
projchp_output(char const *filename,
               int const Nx, int const Ny,
               double const *X, double const *Y,
               double const *U)
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
stationary(char const *filename, struct projchp_method *method)
{
    int Nx = 100, Ny = 100;
    double Lx = 10.0, Ly = 10.0, D = 1.0;

    read_param(filename, &Nx, &Ny, &Lx, &Ly, &D);

    double *U = tdp_vector_new(Nx*Ny);
    double *Uexact = tdp_vector_new(Nx*Ny);
    double *RHS = tdp_vector_new(Nx*Ny);
    double *g = tdp_vector_new(2*Nx); // bottom-top conditions
    double *h = tdp_vector_new(2*Ny); // right-left conditions
    double *X = tdp_vector_new(Nx);
    double *Y = tdp_vector_new(Ny);

    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);

    projchp_grid_init(X, Y, dx, dy, Nx, Ny);

    double B = 2 * D / SQUARE(dx) + 2 * D / SQUARE(dy);
    double Cx = -D / SQUARE(dx);
    double Cy = -D / SQUARE(dy);

    projchp_method_call(method, Nx, Ny, Uexact, X, Y, RHS, g, h, Lx, Ly);

    vector_compute_RHS(Nx, Ny, Cx, Cy, h, g, RHS);
    matrix_5diag_conjugate_gradient(Nx, Ny, B, Cx, Cy, RHS, U);

    projchp_output("numeric.dat", Nx, Ny, X, Y, U);
    projchp_output("exact.dat", Nx, Ny, X, Y, Uexact);

    free(U); free(Uexact);
    free(RHS); free(g); free(h);
    free(X); free(Y);
}

static void
unstationary(char const *filename, struct projchp_method *method)
{
    int Nx = 100, Ny = 100;
    double Lx = 10.0, Ly = 10.0, D = 1.0;
    read_param(filename, &Nx, &Ny, &Lx, &Ly, &D);

    create_directory("sol");
    int const N = Nx * Ny;

    double *U = tdp_vector_new(N);
    double *U0 = tdp_vector_new(N);
    /* double *U4 = tdp_vector_new(N);
       double *U8 = tdp_vector_new(N); */
    double *Uexact = tdp_vector_new(N);
    double *RHS = tdp_vector_new(N);
    double *g = tdp_vector_new(2*Nx); // bottom-top conditions
    double *h = tdp_vector_new(2*Ny); // right-left conditions
    double *X = tdp_vector_new(Nx);
    double *Y = tdp_vector_new(Ny);

    double const Tmax = 10.0;
    int const Nit = 2000;
    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);
    double dt = Tmax / Nit;

    projchp_grid_init(X, Y, dx, dy, Nx, Ny);

    double B = 1.0 / dt + 2.0 * D / SQUARE(dx) + 2.0 * D / SQUARE(dy);
    double Cx = -D / SQUARE(dx);
    double Cy = -D / SQUARE(dy);

    double t = 0.0;
    int const step = 1;
    for (int i = 0; i < Nit; ++i) {
        t += dt;
        method->fu(Nx, Ny, X, Y, RHS, Lx, Ly, t);
        cblas_daxpy(N, 1.0 / dt, U0, 1, RHS, 1);
        vector_compute_RHS(Nx, Ny, Cx, Cy, h, g, RHS);

        matrix_5diag_conjugate_gradient(Nx, Ny, B, Cx, Cy, RHS, U);
        cblas_dcopy(N, U, 1, U0, 1);

        if ((i % step) == 0) {
            char filename[100];
            snprintf(filename, 100, "sol/sol%d.dat", i/step);
            projchp_output(filename, Nx, Ny, X, Y, U);
        }
    }

    printf("Nombre de fichiers: %d\n", Nit/step);
    projchp_output("numeric_unsta.dat", Nx, Ny, X, Y, U);

    free(U); free(Uexact);
    free(U0); //free(U4); free(U8);
    free(RHS); free(g); free(h);
    free(X); free(Y);
}

static struct projchp_method*
projchp_get_method_by_id(unsigned idx)
{
    if (idx >= method_list_length) {
        fprintf(stderr, "No such method!\n");
        exit(EXIT_FAILURE);
    }

    struct projchp_method *l = method_list;
    while (idx--)
        l = l->next;
    return l;
}

static void
solve_equation(struct gengetopt_args_info *opt)
{
    struct projchp_method *m;
    m = projchp_get_method_by_id(opt->function_arg);
    printf("Choosen function: '%s'\n", m->name);
    char const *filename = "param.txt";

    if (m->type == PROJCHP_STATIONARY)
        stationary(filename, m);
    else if (m->type == PROJCHP_UNSTATIONARY)
        unstationary(filename, m);
    else
        projchp_error("invalid method type");
}

static void
handle_opt(struct gengetopt_args_info *opt)
{
    if (opt->list_function_flag) {
        printf("Function list:\n");
        struct projchp_method *list = method_list;
        int i = 0;
        while (list != NULL) {
            printf("%d: %s\n", i, list->name);
            list = list->next;
            ++i;
        }
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char *argv[])
{
    int group_size, rank;
    struct gengetopt_args_info opt;

    cmdline_parser(argc, argv, &opt);
    handle_opt(&opt);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //assert( ((void)"il faut exactement 2 processus MPI", group_size == 2) );

    perf_t p1, p2;
    perf(&p1);
    solve_equation(&opt);
    perf(&p2);
    perf_diff(&p1, &p2);

    uint64_t micro, max_t;
    micro = perf_get_micro(&p2);
    MPI_Reduce(&micro, &max_t, 1, MPI_UNSIGNED_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank)
        fprintf(stderr, "time: %lu.%06lu s\n", max_t/1000000UL, max_t%1000000UL);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
