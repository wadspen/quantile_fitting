#!/bin/bash

#SBATCH --time=08:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=60GB
#SBATCH --constraint=intel
module load gcc
module load r
module load udunits
module load r-rgdal
module load proj
modeul load r gsl

Rscript ./normal_7model_true_mean_wd.R
