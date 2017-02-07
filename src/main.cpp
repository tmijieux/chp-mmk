#include <iostream>

#include "cmdline.h"
#include "proc.hpp"
#include "schwarz_solver.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "func.hpp"
#include "schwarz_printer.hpp"

struct config : public gengetopt_args_info {
    config(int &argc, char **&argv) {
        cmdline_parser(argc, argv, this);
        if (this->list_function_flag) {
            chp::func::print_list();
            throw chp::exception("");
        }
    }
    ~config() { cmdline_parser_free(this); }
};

int main(int argc, char *argv[])
{
    try {
        config opt(argc, argv);
        chp::mpi_proc P;
        chp::timer T;
        chp::schwarz_solver S(P, opt);
        chp::schwarz_printer PR(P, opt);

        T.start();
        S.run(P, PR);
        T.stop();

        T.print(P);

    } catch (chp::exception &e) {
        std::cout << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}
