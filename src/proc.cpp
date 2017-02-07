#include <mpi.h>

#include "proc.hpp"

using namespace chp;

mpi_proc::mpi_proc()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &m_group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
}

mpi_proc::~mpi_proc()
{
    MPI_Finalize();
}
