#!/bin/bash
#PBS -N submit_fa2domain
#PBS -l nodes=1:ppn=1 (FA: originally 8)
## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Joe Burke
OUTPATH=$PBS_O_WORKDIR
cd ${OUTPATH}
## USER INPUTS
INFILE=$1
TYPE=$2
PHYLO=$3
DIR=$4
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}
echo "Splitting queries into domains, please wait" > status.txt
export INTERPRO=/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:/var/interproscan/data:$INTERPRO
export SIGNALP=/var/interproscan/bin/signalp/4.1
source /etc/profile.d/modules.sh
module purge
module load iprscan
module load R
# FA: disabled b/c docker's handling perms
# setfacl -R -m group:shiny:r-x ${DIR}
sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${INFILE} ${PREFIX} ${DIR}
touch ${DIR}/${PREFIX}.domains.fa
Rscript --vanilla /data/research/jravilab/molevol_scripts/upstream_scripts/fa2domain.R $INFILE ${DIR}/${PREFIX}.iprscan.tsv ${DIR}/${PREFIX}.domains.fa ${DIR} ${PHYLO} ${TYPE}