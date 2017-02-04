#include <stdio.h>
#include <stdlib.h>

#include "cmdline.h"
#include "proc.h"
#include "schwarz_solver.h"
#include "schwarz_printer.h"
#include "timer.h"
#include "func.h"

static void handle_opt(struct gengetopt_args_info *opt)
{
    if (opt->list_function_flag) {
        chp_print_func_list();
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char *argv[])
{
    struct gengetopt_args_info opt;
    chp_proc P;
    chp_schwarz_solver S;
    chp_timer T;
    chp_schwarz_printer PR;

    cmdline_parser(argc, argv, &opt);
    handle_opt(&opt);

    chp_proc_init(&P);
    chp_schwarz_printer_init(&PR, &P, opt.verbose_flag, opt.file_flag);
    chp_schwarz_solver_init(&S, &P, &opt);

    chp_timer_start(&T);
    chp_schwarz_solver_run(&S, &P, &PR);
    chp_timer_stop(&T);

    chp_timer_print(&T, &P, opt.recouvr_arg);

    chp_schwarz_solver_free(&S);
    cmdline_parser_free(&opt);
    chp_proc_fini(&P);
    return EXIT_SUCCESS;
}
