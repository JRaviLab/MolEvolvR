#!/bin/bash

#PBS -l nodes=1:ppn=1          # number of nodes requested (FA: originally 10)
#PBS -m abe                     # email notifications for job
#PBS -M=sosinsk7@msu.edu        # user email; RESET

# test segfault error in ipr2da

module load R

Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R $1 $2 $3 $4

printf "\nFinished.\n"
