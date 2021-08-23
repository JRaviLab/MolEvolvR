#!/bin/bash
#PBS -N wrapper_fa2domain
#PBS -l nodes=1 :ppn 8
## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Joe Burke

## USER INPUTS
echo $INFILE>> log.new
echo $SCRIPT>> log.new
echo $OTHERARGS >> log.new
DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}
export INTERPRO=/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:$INTERPRO
source /etc/profile.d/modules.sh
module purge
module load iprscan
module load R
setfacl -R -m group:shiny:r-x ${DIR}
sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${INFILE} ${PREFIX} ${DIR}
touch ${DIR}/${PREFIX}.domains.fa
Rscript --vanilla /data/research/jravilab/molevol_scripts/upstream_scripts/fa2domain.R $INFILE ${DIR}/${PREFIX}.iprscan.tsv ${DIR}/${PREFIX}.domains.fa

sh ${SCRIPT} ${DIR}/${PREFIX}.domains.fa ${OTHERARGS}
