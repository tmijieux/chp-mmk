#include <iostream>
#include <sstream>
#include "schwarz_printer.hpp"

using namespace chp;


static void output(const char *filename, int const Nx, int const Ny,
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

schwarz_printer::schwarz_printer(const proc& p, bool console_output, bool file_output):
    m_proc(p),
    m_console_output(console_output),
    m_file_output(file_output)
{
}

schwarz_printer::schwarz_printer(const proc &p, struct gengetopt_args_info const& opt):
    schwarz_printer(p, opt.verbose_flag, opt.file_flag)
{
}

void
schwarz_printer::
print_unsta_step(int step, int schwarz_step, const equation& eq) const
{
    if (m_console_output && !m_proc.rank())
        std::cout << "# time_step("<< step << "): "
                  << "schwarz_step = "<< schwarz_step
                  << std::endl;

    if (m_file_output && (step % print_freq) == 0) {
        std::stringstream ss;
        ss << "sol/sol" << step/print_freq + 1 <<".dat." <<  m_proc.rank();
        output(ss.str().c_str(), eq.Nx, eq.Ny, eq.X, eq.Y, eq.U1);
    }
}

void
schwarz_printer::
print_unsta_fin(int total_step, int total_schwarz_step) const
{
    if (m_console_output && !m_proc.rank())
        std::cout << "# total_schwarz_step = " << total_schwarz_step << std::endl
                  << "# total_step = "<< total_step << std::endl;
}

void
schwarz_printer::
print_sta_step(int schwarz_step, int step) const
{
    if (m_console_output && !m_proc.rank())
        std::cout << "schwarz_step(" << schwarz_step << "): "
                  << "solver_step = " << step << std::endl;
}

void
schwarz_printer::
print_sta_fin(int schwarz_step, int total_step, const equation &eq) const
{
    if (m_console_output && !m_proc.rank())
        std::cout << "# schwarz_step = "<< schwarz_step << std::endl
                  << "# total_step = " << total_step << std::endl;

    if (m_file_output) {
        std::stringstream ss;
        ss << "numeric.dat." << m_proc.rank();
        output(ss.str().c_str(), eq.Nx, eq.Ny, eq.X, eq.Y, eq.U1);
    }
}
