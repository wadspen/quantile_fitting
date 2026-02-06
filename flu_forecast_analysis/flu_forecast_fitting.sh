#!/bin/bash

#SBATCH --time=5-10:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --constraint=intel
#module purge
#module load r/4.4.1 gcc udunits proj gsl r-rgdal
module purge
module unuse /opt/rit/el9/20250815/modules/lmod/linux-rhel9-x86_64/Core/
module load r/4.4.1 gcc udunits proj gsl r-rgdal


Rscript ./final_analysis.R "$1" "$2" "$3" "$4"

