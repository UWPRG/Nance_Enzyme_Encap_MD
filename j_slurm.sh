#!/bin/bash 

## Job Name 
#SBATCH --job-name=XXX_wat_system

## Allocation Definition
#SBATCH --account=pfaendtner
#SBATCH --partition=ckpt

## Resources 
## Nodes 
#SBATCH --nodes=XX       

## Tasks per node (Slurm assumes you want to run 28 tasks, remove 2x # and adjust parameter if needed)
###SBATCH --ntasks-per-node=16 

## Walltime (ten minutes) 
#SBATCH --time=XX 

# E-mail Notification, see man sbatch for options
#SBATCH --mail-type=NONE

## Memory per node 
#SBATCH --mem=44G 

## Specify the working directory for this job 
#SBATCH --workdir=XXX

# Loading 2016 gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

## loading plumed
export PATH="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/src/lib/:$PATH"
export LIBRARY_PATH="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/src/lib/:$LIBRARY_PATH"
export LD_LIBRARY_PATH="/gscratch/pfaendtner/cdf6gc/codes/libmatheval/lib/.libs/:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/src/lib/:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/src/lib/:$DYLD_LIBRARY_PATH"
export PLUMED_KERNEL="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/src/lib/libplumedKernel.so"
export PLUMED_VIMPATH="/suppscr/pfaendtner/arushi3/pbmetad-ip/plumed2/vim"

##module load icc_15.0.3-impi_5.0.3

mpirun -np XX gmx_mpi mdrun -cpi restart -cpo restart -append -cpt 1.0 &> log.txt

#./YYY.sh

