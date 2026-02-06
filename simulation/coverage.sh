#!/bin/bash

#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=80GB
#SBATCH --constraint=intel

module purge
module unuse /opt/rit/el9/20250815/modules/lmod/linux-rhel9-x86_64/Core/
module load r/4.4.1 gcc udunits proj gsl r-rgdal


#Rscript ./variational_7model_uwd.R "$1" $2 $3
Rscript ./normal_7model_true_mean.R "$1" "$2" "$3"
#Rscript ./semipar_cover.R "$1" "$2" "$3"
