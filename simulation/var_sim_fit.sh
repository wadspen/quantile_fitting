#!/bin/bash

#SBATCH --time=5-23:59:00
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

#Rscript ./variational_7model_uwd.R "$1" $2 $3
#Rscript ./normal_7model_wd.R "$1" "$2" "$3"
#Rscript ./exponential_7model_wd.R "$1" "$2" "$3"
Rscript ./tukey_2model_time.R "$1" "$2"
