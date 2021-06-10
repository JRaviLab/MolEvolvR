#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Joe Burke

## USER INPUTS
INFILE=$1
SCRIPT=$2
OTHERARGS=$3
DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}

sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh $INFILE ${PREFIX} $DIR

Rscript --vanilla fa2domain.R $INFILE ${PREFIX}.iprscan domains.fa

qsub ${SCRIPT} domains.fa ${OTHERARGS}
