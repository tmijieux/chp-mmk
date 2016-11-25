#!/usr/bin/env bash
#SBATCH --job-name=mijieux0
#SBATCH --output=out.0
#SBATCH --error=err.0
#SBATCH -p mistral
#SBATCH --time=02:00:00
#  # SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 10

# 4 noeud / 1 proc mpi par noeud

WORKDIR=${WORKDIR:-${HOME}/chp-mmk}
MPIEXEC=mpirun

cd ${WORKDIR}
. ./.module.load

mpd & 

do_job() {
    size=$1
    ${MPIEXEC} -n $size ./projCHP
    rm $file
}

for i in $(seq 1 20); do
    do_job $i
done

