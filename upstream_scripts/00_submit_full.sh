#!/bin/bash
#PBS -N submit_full
#PBS -l nodes=1 :ppn 2
## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## USER INPUTS
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

if [ $WBLAST = T ]; then
   mv "wblast.tsv" "${PREFIX}.wblast.tsv"
   find $PWD -type f -name "${PREFIX}.wblast.tsv" > input.txt
   INPATHS=input.txt
   echo "qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F '$INPATHS $WBLAST'" >> cmd.txt
   ID=`qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "$INPATHS $WBLAST"`
fi

if [ $WBLAST = F ]; then

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
      printf "No. of seqs provided: $NFASTA\nSo, we are going to proceed to the analysis.\n"
      find $PWD -type f -name $BASE > input.txt
   fi

   INPATHS=input.txt
   echo "0/${NFASTA} analyses completed" > status.txt 
   ID=`qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "$INPATHS $WBLAST" -t 1-$NFASTA`

fi

setfacl -R -m group:shiny:r-x ${DIR}
#NUM_ID=`echo ${ID} | grep -Eo "^[[:digit:]]+"`
#run_start=0
#qstat | grep ${NUM_ID}
#while [ $? -ne 1 ] && [ $run_start -lt 86400 ];
#do
#  sleep 3600
#  ((run_start+=3600))
#  qstat | grep ${NUM_ID}
#done
#if [ $run_start -gt 86400 ];
#   then
#   qdel ${ID}
#   touch "kill.txt"
#   exit 1
#fi
#exit 0
