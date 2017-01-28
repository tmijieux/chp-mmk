#!/usr/bin/env bash
#SBATCH --job-name=mijieux0
#SBATCH --output=out.0
#SBATCH --error=err.0
#SBATCH -p mistral
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 20

# 4 noeud / 1 proc mpi par noeud

WORKDIR=${WORKDIR:-${HOME}/chp-mmk}
MPIEXEC=mpiexec

cd ${WORKDIR}
. ./.module.load

export MKL_NUM_THREADS=20
${MPIEXEC} -np 4 ./projCHP -X 100 -Y 100 -f 0 -x 1.0 -y 1.0 -R 5

