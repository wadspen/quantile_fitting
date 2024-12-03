#!/bin/bash

#SBATCH --time=7-10:59:00
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
modeul load r gsl

Rscript ./flu_forecast_fitting.R
