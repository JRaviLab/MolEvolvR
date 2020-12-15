#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## Input type: fasta file
## Output type: N individual fasta files, text file containing paths to individual fasta files

## USER INPUTS
INFILE=$1
DB=refseq
NHITS=5000
EVAL=1e-5

dir=$(dirname $INFILE)
base=$(basename $INFILE)
cd ${dir}

printf "START_DT\tSTOP_DT\tquery\tdblast_duration\tdblast_clnup_duration\tacc2fa_duration\tbclust_duration\tc2t_duration\tipr_duration}\tipr2da_duration\trps_duration\trps2da_duration\tduration\n" >> ${dir}/logfile.txt

## load modules ##
module load R
module load iprscan

## query protein runs
query_run_start=$SECONDS
sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${INFILE} $base $dir
Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${base}.iprscan.tsv $base
Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${base}.iprscan.tsv $base
query_run_duration=$(( $SECONDS - $query_run_start ))

printf "\nQuery run duration: $query_run_duration\n"
#printf "$query_run_duration"  >> ${dir}/logfile.txt
splitfa_start=$SECONDS

NFASTA=$(grep -o "^>" $INFILE | wc -l)

if [ $NFASTA -gt 1 ]
then
   printf "No. of seqs provided: $NFASTA>1\nSo, we are going to split it up for you prior to the analysis.\n"
     # https://unix.stackexchange.com/questions/15662/splitting-text-files-based-on-a-regular-expression
     awk -F "( )|(>)" '/^>/{x=""$2".faa";}{print >x;}' $INFILE
   #sh /data/research/jravilab/molevol_scripts/upstream_scripts/split_fasta.sh ${INFILE}
   find $PWD -type f -name "*.faa" > input.txt
fi

if [ $NFASTA = 1 ] 
then
   printf "No. of seqs provided: $NFASTA=1\nSo, we are going to proceed to the analysis.\n"
   find $PWD -type f -name "$base" > input.txt
fi

splitfa_duration=$(( $SECONDS - $splitfa_start ))
printf "\nSplit FASTA file duration: $splitfa_duration\n"
#printf "$splitfa_duration"  >> ${dir}/logfile.txt

INPATHS=input.txt
#cat $INPATHS
qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F $INPATHS

chmod ug+rw -R * 
