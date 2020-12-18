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

DIR=$(dirname $INFILE)
BASE=$(basename $INFILE)
PREFIX=$(echo "${BASE%%.*}")
cd ${DIR}

printf "START_DT\tSTOP_DT\tquery\tdblast_duration\tacc2fa_duration\tacc2taxlin_dur\tdblast_clnup_duration\tbclust_duration\tc2t_duration\tipr_duration}\tipr2da_duration\trps_duration\trps2da_duration\tduration\n" >> ${DIR}/logfile.txt

## LOAD MODULES ##
module load R
module load edirect
module load iprscan

## QUERY PROT ONLY RUNS ##
query_run_start=$SECONDS

## if wblast == T, take 1st protein and run ipscan analysis

if [ $WBLAST = T ]; then
   sed 's/,/\t/g' ${INFILE} >> ${OUTDIR}/${PREFIX}.wblast.tsv

   QUERY=$(cat ${OUTDIR}/${PREFIX}.wblast.tsv | head -n 1 | awk '{print $1}')
   printf "$QUERY"
   efetch -db protein -format fasta -id "$QUERY" >> ${OUTDIR}/${QUERY}.fa

   ## IPRSCAN ##
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${OUTDIR}/${QUERY}.fa ${QUERY} $DIR
   # grab accnum from iprscan file, append to all_accnum file
   cat ${DIR}/${BASE}.iprscan.tsv | awk -F "\t" '{ print $1 }' | sort -u >> ${DIR}/${QUERY}.all_accnums.txt
   ## ACC2TAX ##
   printf "\nACC2TAX\n"
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2tax.sh ${DIR}/${QUERY}.all_accnums.txt ${QUERY} ${DIR}
   ## IPR2LIN
   printf "\nIPR2LIN\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${DIR}/${QUERY}.iprscan.tsv ${DIR}/${QUERY}.acc2info.tsv $QUERY
   ## IPR2DA ##
   printf "\nIPR2DA\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${DIR}/${QUERY}.iprscan.tsv $QUERY

   query_run_duration=$(( $SECONDS - $query_run_start ))

   printf "\nQuery run duration: $query_run_duration\n"
   splitfa_start=$SECONDS

   find $PWD -type f -name "$QUERY" > input.txt
fi

if [ $WBLAST = F ]; then

   ## IPRSCAN ##
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${INFILE} $PREFIX $DIR
   # grab accnum from iprscan file, append to all_accnum file
   cat ${DIR}/${PREFIX}.iprscan.tsv | awk -F "\t" '{ print $1 }' | sort -u >> ${DIR}/${PREFIX}.all_accnums.txt
   ## ACC2TAX ##
   printf "\nACC2TAX\n"
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2tax.sh ${DIR}/${PREFIX}.all_accnums.txt ${PREFIX} ${DIR}
   ## IPR2LIN
   printf "\nIPR2LIN\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${DIR}/${PREFIX}.iprscan.tsv ${DIR}/${PREFIX}.acc2info.tsv $PREFIX
   ## IPR2DA ##
   printf "\nIPR2DA\n"
   Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${DIR}/${PREFIX}.iprscan.tsv $PREFIX

   query_run_duration=$(( $SECONDS - $query_run_start ))

   printf "\nQuery run duration: $query_run_duration\n"
   splitfa_start=$SECONDS

   ## COUNT NUMBER OF FASTA SEQUENCES ##
   NFASTA=$(grep -o "^>" $INFILE | wc -l)

   if [ $NFASTA -gt 1 ]
   then
      printf "No. of seqs provided: $NFASTA>1\nSo, we are going to split it up for you prior to the analysis.\n"
      # https://unix.stackexchange.com/questions/15662/splitting-text-files-BASEd-on-a-regular-expression
      awk -F "( )|(>)" '/^>/{x=""$2".faa";}{print >x;}' $INFILE
      find $PWD -type f -name "*.faa" > input.txt
   fi

   if [ $NFASTA = 1 ] 
   then
      printf "No. of seqs provided: $NFASTA\nSo, we are going to proceed to the analysis.\n"
      find $PWD -type f -name "$BASE" > input.txt
   fi

   splitfa_duration=$(( $SECONDS - $splitfa_start ))
   printf "\nSplit FASTA file duration: $splitfa_duration\n"
fi 

INPATHS=input.txt
#qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "$INPATHS $WBLAST"

