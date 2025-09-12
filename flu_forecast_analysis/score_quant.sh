#!/bin/bash

#SBATCH --time=3-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel
module purge
module load r/4.4.1 gcc udunits proj gsl r-rgdal


Rscript ./score_quantile_fits.R "$1" 

