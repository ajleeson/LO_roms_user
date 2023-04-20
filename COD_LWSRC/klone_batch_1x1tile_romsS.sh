#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=1

## Walltime 
#SBATCH --time=00:20:00

## Memory per node
#SBATCH --mem=128G

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort/
export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

RUN_DIR=/mmfs1/gscratch/macc/auroral/LO_roms_user/COD_LWSRC
# mpirun -np 1 $RUN_DIR/romsS $RUN_DIR/roms_cod_LwSrc_1x1tile_romsS.in > $RUN_DIR/roms_log_1x1_romsS.txt
$RUN_DIR/romsS $RUN_DIR/roms_cod_LwSrc_1x1tile_romsS.in > $RUN_DIR/roms_log_1x1_romsS.txt

