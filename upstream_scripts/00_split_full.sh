#!/bin/bash

## USER INPUTS
INFILE=$1
DB=refseq
NHITS=5000
EVAL=1e-5

dir=$(dirname $INFILE)
base=$(basename $INFILE)
cd ${dir}


NFASTA=$(grep -o "^>" input.txt | wc -l)


if [ $NFASTA -gt 1 ]
then
   printf "No. of seqs provided: $NFASTA>1\nSo, we are going to split it up for you prior to the analysis.\n"
   sh /data/research/jravilab/molevol_scripts/upstream_scripts/split_fasta_jr.sh ${INFILE}
   find $PWD -type f -name "*.faa" > input.txt
fi

if [ $NFASTA = 1 ] 
then
   printf "No. of seqs provided: $NFASTA=1\nSo, we are going to proceed to the analysis.\n"
   find $PWD -type f -name "$base" > input.txt
fi

INPATHS=input.txt
cat $INPATHS
qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper.sb -F $INPATHS
