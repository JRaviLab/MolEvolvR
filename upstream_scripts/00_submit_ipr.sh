#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## Input type: web or command-line IPRSCAN runs | expects iprscan.tsv files
## Output type: IPR file with DomArch added + protein metadata such as TaxID, Species, Lineage, Name & such

## USER INPUTS

## USAGE
## ./00_submit_ipr.sh accnum_iprscan.tsv

DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
cd ${DIR}

find $PWD -type f -name "$BASE" > input.txt

INPATHS=input.txt
echo "0/1 analyses completed" > status.txt
qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb -F $INPATHS

# FA: disabled b/c docker's handling perms
# setfacl -R -m group:shiny:r-x ${DIR}
