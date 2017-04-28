#!/bin/bash
#
#SBATCH --account=plantanalytics
#SBATCH --time=1:00
#SBATCH --nodes=3
#SBATCH --output=outParallel.out

module load intel
module load R/3.3.1

R < parallelMtMoran.R --no-save 
