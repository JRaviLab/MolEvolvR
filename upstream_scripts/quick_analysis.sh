#!/bin/bash

## USER INPUTS
INFILE=$1
FILE_TYPE=$2

## USAGE ##
# quick analysis on query proteins/input files
# takes in web blast or FASTA file

DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}

## LOAD MODULES ##
module load R
module load edirect
module load iprscan

## QUERY PROT ONLY RUNS ##
query_run_start=$SECONDS

mkdir ${PREFIX}_quick_out
cp ${INFILE} ${PREFIX}_quick_out
cd ${PREFIX}_quick_out

if [ $FILE_TYPE = WBLAST ]; then
   
   sed 's/,/\t/g' ${INFILE} > ${DIR}/${PREFIX}.wblast.tsv
   QUERY=$(cat ${DIR}/${PREFIX}.wblast.tsv | head -n 1 | awk '{ print $1 }')
   printf "$QUERY" > ${DIR}/${QUERY}.accnum.txt
   efetch -db protein -format fasta -id $QUERY > ${DIR}/${QUERY}.faa
   
   ## IPRSCAN ##
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${DIR}/${QUERY}.faa ${QUERY} $DIR
      
   ## ACC2INFO ##
   printf "\nACC2INFO\n"
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh ${DIR}/${QUERY}.accnum.txt ${QUERY} ${DIR}
   
   ## IPR2LIN
   printf "\nIPR2LIN\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${DIR}/${QUERY}.iprscan.tsv ${DIR}/${QUERY}.acc2info.tsv $QUERY

   ## IPR2DA ##
   printf "\nIPR2DA\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${DIR}/${QUERY}.iprscan_cln.tsv $QUERY

   query_run_duration=$(( $SECONDS - $query_run_start ))
   printf "\nQuery run duration: $query_run_duration\n"

fi


if [ $FILE_TYPE = FASTA ]; then

   ## IPRSCAN ##
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${INFILE} $PREFIX $DIR         
   
   # grab accnum from iprscan file, append to all_accnum file
   cat ${DIR}/${PREFIX}.iprscan.tsv | awk -F "\t" '{ print $1 }' | sort -u >> ${DIR}/${PREFIX}.accnums.txt
   
   ## ACC2INFO ##
   printf "\nACC2INFO\n"
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh ${DIR}/${PREFIX}.accnums.txt ${PREFIX} ${DIR}
   
   ## IPR2LIN
   printf "\nIPR2LIN\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${DIR}/${PREFIX}.iprscan.tsv ${DIR}/${PREFIX}.acc2info.tsv $PREFIX
   
   ## IPR2DA ##
   printf "\nIPR2DA\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${DIR}/${PREFIX}.iprscan.tsv ${DIR}/${PREFIX}.acc2info.tsv $PREFIX
   
   query_run_duration=$(( $SECONDS - $query_run_start ))
   printf "\nQuery run duration: $query_run_duration\n"
   
fi

setfacl -R -m group:shiny:r-x ${DIR}
