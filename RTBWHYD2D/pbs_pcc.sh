#!/bin/bash
#PBS -N BlastwaveRT
#PBS -q long
#PBS -m n
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -j oe

module load intel/2024

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd $PBS_O_WORKDIR

./Simulation.x

