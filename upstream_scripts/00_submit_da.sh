#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## Input type: fasta file
## Output type: text file containing paths to individual fasta files

## USER INPUTS
INFILE=$1

## USAGE
## ./00_submit_da.sh multifasta.fa

DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
cd ${DIR}

printf "START_DT\tSTOP_DT\tquery\tacc2info\tblast_clust\tclust2table\tiprscan\tipr2da\tduration\n" >> ${DIR}/logfile.txt

# grep "|" handles files that are not in ncbi format
# split each word in the header by "|" and use the second element as the accNum
grep "|" $INFILE
if [ $? = 0 ]
then
awk -F '|' '/^>/{x=""$2".faa";}{print >x;}' $INFILE
  find $PWD -type f -name "*.faa" > input.txt
else
find $PWD -type f -name $BASE > input.txt
fi



INPATHS=input.txt
NFASTA=$(grep -o "^>" $INFILE | wc -l)
echo "0/${NFASTA} analyses completed" > status.txt 
qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_da.sb -F $INPATHS

setfacl -R -m group:shiny:r-x ${DIR}
