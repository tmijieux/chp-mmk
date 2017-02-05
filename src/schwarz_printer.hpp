#ifndef CHP_SCHWARZ_PRINTER_H
#define CHP_SCHWARZ_PRINTER_H

#include "proc.hpp"
#include "equation.hpp"
#include "cmdline.h"


namespace chp {

class schwarz_printer{
    static const int print_freq = 100;
private:
    const proc& m_proc;
    bool m_console_output;
    bool m_file_output;


public:

    schwarz_printer(const proc& p, bool console_output = false, bool file_output = false);
    schwarz_printer(const proc& p, struct gengetopt_args_info const &);

    void print_unsta_step(int step, int schwarz_step, const equation& eq) const;
    void print_unsta_fin(int total_step, int total_schwarz_step) const;

    void print_sta_step(int schwarz_step, int step) const;
    void print_sta_fin(int schwarz_step, int total_step, const equation &eq) const;
};


};
#endif // CHP_SCHWARZ_PRINTER_H
