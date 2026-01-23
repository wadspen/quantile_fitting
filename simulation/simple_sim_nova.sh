#!/bin/bash

#SBATCH --time=7:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel

#module purge
#module load zlib-ng/2.1.6-nf73pqi curl/8.11.1-xzf3njp proj/9.4.1-3qfgrvd gcc udunits gsl
#module load r/4.2.2-py310-ly4mhww #r/4.4.1 #r-gdalutils/2.0.3.2-py310-r42-openmpi4-7loeoid
#module load r-rgdal
module purge
module unuse /opt/rit/el9/20250815/modules/lmod/linux-rhel9-x86_64/Core/
module load r/4.4.1 gcc udunits proj gsl r-rgdal



Rscript simple_sim_cover.R "$1" "$2" "$3" "$4" "$5"
