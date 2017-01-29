#include <mpi.h>

#include "proc.hpp"

using namespace chp;

proc::proc()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &m_group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
}

proc::~proc()
{
    MPI_Finalize();
}
