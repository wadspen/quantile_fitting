#!/bin/bash

#SBATCH --time=5-10:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel
module purge
module load r/4.4.1 gcc udunits proj gsl r-rgdal



Rscript ./flu_forecast_fitting_ord.R "$1"

