#!/bin/bash
#PBS -N submit_fa2domain
#PBS -l nodes=1 :ppn 8
## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Joe Burke

## USER INPUTS
echo $INFILE
echo $SCRIPT
echo $OTHERARGS
DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}
module load iprscan
module load R

sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh $INFILE ${PREFIX} $DIR

Rscript --vanilla /data/research/jravilab/molevol_scripts/upstream_scripts/fa2domain.R $INFILE ${PREFIX}.iprscan ${DIR}/domains.fa

sh ${SCRIPT} domains.fa ${OTHERARGS}
