#include "schwarz_printer.h"

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

void chp_schwarz_printer_stationary(
    chp_schwarz_printer *pr, int schwarz_step, int total_step, const chp_equation *eq)
{
    if (pr->verbose_output) {
        if (!pr->p->rank)
            printf("# schwarz_step = %d\n# total_step = %d\n", schwarz_step, total_step);

    }
    if (pr->file_output) {
        char filename[100];
        snprintf(filename, 100, "numeric.dat.%d", p->rank);
        chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
    }
}

void chp_schwarz_printer_unstationary(
    chp_schwarz_printer *pr, int step, const chp_equation *eq)
{
    if (pr->file_output) {
        if ((step % pr->print_freq) == 0) {
            char filename[100];
            snprintf(filename, 100,
                     "sol/sol%d.dat.%d", i/print_step + 1, pr->p->rank);
            chp_output(filename, eq->Nx, eq->Ny, eq->X, eq->Y, eq->U1);
        }
    }

}

void chp_schwarz_printer_unstationary_final(
    chp_schwarz_printer *pr, int schwarz_step, int total_step)
{
    if (pr->verbose_output && !pr->p->rank) {
        printf("# schwarz_step = %d\n# total_step = %d\n", schwarz_step, total_step);
}


void chp_schwarz_printer_init(
    chp_schwarz_printer*pr, const chp_proc*p, bool verbose, bool file)
{
    memset(pr, 0 sizeof*pr);
    pr->verbose_output = verbose;
    pr->file_output = file;
    pr->p = p;
    pr->print_freq = 100;
}


