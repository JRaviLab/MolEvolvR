#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## USER INPUTS
INFILE=$1
DB=refseq
NHITS=5000
EVAL=1e-5

## USAGE ##
# Full analysis | input | fasta file
# sh /path/to/00_submit_full.sh /path/to/multifasta.fa F
# Web-BLAST analysis | input | blastp file as csv
# sh /path/to/00_submit_full.sh /path/to/web_blast.csv T

DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}

## append column names to logfile, data gets appended in wrapper scripts
printf "START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration\n" >> ${DIR}/logfile.tsv

## QUERY PROT ONLY RUNS ##
query_run_start=$SECONDS

## COUNT NUMBER OF FASTA SEQUENCES ##
NFASTA=$(grep -o "^>" $INFILE | wc -l)

if [ $NFASTA -gt 1 ]
then
  printf "No. of seqs provided: $NFASTA>1\nSo, we are going to split it up for you prior to the analysis.\n"
  # https://unix.stackexchange.com/questions/15662/splitting-text-files-BASEd-on-a-regular-expression
  # grep "|" handles files that are not in ncbi format
  # split each word in the header by "|" and use the second element as the accNum
  grep "|" $INFILE
  if [ $? = 0 ]
  then
awk -F '|' '/^>/{x=""$2".faa";}{print >x;}' $INFILE
    find $PWD -type f -name "*.faa" > input.txt
else
  awk -F "( )|(>)" '/^>/{x=""$2".faa";}{print >x;}' $INFILE
  find $PWD -type f -name "*.faa" > input.txt
fi
fi

if [ $NFASTA = 1 ]
then
  grep "|" $INFILE
  if [ $? = 0 ]
  then
awk -F '|' '/^>/{x=""$2".faa";}{print >x;}' $INFILE
  fi
  find $PWD -type f -name "*.faa" > input.txt
fi

INPATHS=input.txt

qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F "$INPATHS" -t 1-$NFASTA

setfacl -R -m group:shiny:r-x ${DIR}
