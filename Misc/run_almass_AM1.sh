#!/bin/bash
#SBATCH --account almass_genes_ana
#SBATCH -c 16
#SBATCH --mem 32g
#SBATCH --time 50:00:00
singularity exec singularity.sif python3 my_running_ALMaSS_script_AM.py
