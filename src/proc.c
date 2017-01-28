#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "proc.h"

void chp_proc_init(chp_proc *P)
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &P->group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &P->rank);
}

void chp_proc_fini(chp_proc *P)
{
    memset(P, 0, sizeof*P);
    MPI_Finalize();
}
