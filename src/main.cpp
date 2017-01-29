#include <iostream>

#include "cmdline.h"
#include "proc.hpp"
#include "schwarz_solver.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "func.hpp"

using namespace chp;


struct config : public gengetopt_args_info {
    config(int &argc, char **&argv) {
        cmdline_parser(argc, argv, this);
        if (this->list_function_flag) {
            func::print_list();
            throw exception("");
        }
    }
    ~config() { cmdline_parser_free(this); }
};

int main(int argc, char *argv[])
{
    try {
        config opt(argc, argv);
        proc P;
        timer T;
        schwarz_solver S(P, &opt);

        T.start();
        S.run(P);
        T.stop();

        T.print(P);

    } catch (exception &e) {
        std::cout << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}
