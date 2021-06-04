#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## USER INPUTS
INFILE=$1
DB=refseq
NHITS=5000
EVAL=1e-5
WBLAST=$2

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

if [ $WBLAST = T ]; then
   sed 's/,/\t/g' ${INFILE} > ${DIR}/${PREFIX}.wblast.tsv
   find $PWD -type f -name "${PREFIX}.wblast.tsv" > input.txt

   INPATHS=input.txt

   qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "$INPATHS $WBLAST"
fi

if [ $WBLAST = F ]; then

   ## COUNT NUMBER OF FASTA SEQUENCES ##
   NFASTA=$(grep -o "^>" $INFILE | wc -l)

   if [ $NFASTA -gt 1 ]
   then
      printf "No. of seqs provided: $NFASTA>1\nSo, we are going to split it up for you prior to the analysis.\n"
      # https://unix.stackexchange.com/questions/15662/splitting-text-files-BASEd-on-a-regular-expression
      grep "|" $INFILE
      if [ $? ]
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
      printf "No. of seqs provided: $NFASTA\nSo, we are going to proceed to the analysis.\n"
      grep "|" $INFILE
      if [ $? ]
      then
	awk -F '|' '/^>/{x=""$2".faa";}{print >x;}' $INFILE
      fi
      find $PWD -type f -name "*.faa" > input.txt
   fi

   INPATHS=input.txt

   qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "$INPATHS $WBLAST" -t 1-$NFASTA

fi

setfacl -R -m group:shiny:r-x ${DIR}
