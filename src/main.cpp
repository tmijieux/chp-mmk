#include <iostream>

#include "cmdline.h"
#include "proc.hpp"
#include "schwarz_solver.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "func.hpp"

using namespace chp;

static void handle_opt(struct gengetopt_args_info *opt)
{
    if (opt->list_function_flag) {
        func::print_list();
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char *argv[])
{
    try {
        struct gengetopt_args_info opt;
        proc P;
        timer T;
        
        cmdline_parser(argc, argv, &opt);
        handle_opt(&opt);
        
        schwarz_solver S(P, &opt);

        T.start();
        S.run(P);
        T.stop();
        
        T.print(P);
        
        cmdline_parser_free(&opt);
    } catch (exception &e) {
        std::cout << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}
