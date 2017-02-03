#ifndef CHP_SCHWARZ_PRINTER_H
#define CHP_SCHWARZ_PRINTER_H

#include <stdbool.h>

typedef struct chp_schwarz_printer_ chp_schwarz_printer;

#include "proc.h"
#include "equation.h"

struct chp_schwarz_printer_ {
    const chp_proc *p;
    int print_freq;
    bool verbose_output;
    bool file_output;
};


void chp_schwarz_printer_init(chp_schwarz_printer*pr, const chp_proc*p,
                              bool verbose, bool file);

void chp_schwarz_printer_stationary(
    chp_schwarz_printer *pr, int schwarz_step, int total_step, const chp_equation *eq);
void chp_schwarz_printer_unstationary(
    chp_schwarz_printer *pr, int step, const chp_equation *eq);
void chp_schwarz_printer_unstationary_final(
    chp_schwarz_printer *pr, int schwarz_step, int total_step);

#endif // CHP_SCHWARZ_PRINTER_H
