#!/bin/bash
INFILE=$1
SCRIPT=$2
OTHERARGS=$3

qsub -v SCRIPT=${SCRIPT},INFILE=${INFILE},OTHERARGS=${OTHERARGS} /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_fa2domain.sh
