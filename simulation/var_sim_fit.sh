#!/bin/bash

#SBATCH --time=5-23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel

module purge
module load r/4.4.1 gcc udunits proj gsl r-rgdal

#Rscript ./variational_7model_uwd.R "$1" $2 $3
#Rscript ./normal_7model_wd.R "$1" "$2" "$3"
#Rscript ./exponential_7model_wd.R "$1" "$2" "$3"
Rscript ./var_num_comps_uwd.R "$1" "$2" "$3"
