#!/bin/bash

#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel

rm -rf ./sim_draws
