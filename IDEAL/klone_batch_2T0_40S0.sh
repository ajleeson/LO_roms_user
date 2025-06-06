#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=4

## Walltime 
#SBATCH --time=00:20:00

## Memory per node
#SBATCH --mem=128G

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort/
export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

RUN_DIR=/mmfs1/gscratch/macc/auroral/LO_roms_user/IDEAL
mpirun -np 1 $RUN_DIR/romsM $RUN_DIR/roms_2T0_40S0.in > $RUN_DIR/2T0_40S0_roms_log.txt

