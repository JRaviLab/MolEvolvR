#!/bin/bash
#PBS -N submit_fa2domain
#PBS -l nodes=1 :ppn 1
qsub -v SCRIPT=${SCRIPT},INFILE=${INFILE},OTHERARGS=${OTHERARGS} /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_fa2domain.sh
