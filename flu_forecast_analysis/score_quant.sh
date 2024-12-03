#!/bin/bash

#SBATCH --time=2-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel
module load gcc
module load r
module load udunits
module load r-rgdal
module load proj
module load r gsl


Rscript ./score_quantile_fits.R "$1" 

