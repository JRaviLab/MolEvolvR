#!/bin/bash

# created by Lauren Sosinski 2019.01.21
# last edit 2020.06.25

################
## THE SCRIPT ##
################

## FILEPATHS: I/O
INFILE=$1
suffix=$2
outdir=$3

OUTFILE=$(printf "${outdir}/${suffix}.rps.out")

printf "\n###################################"
printf "\n Now running RPSBLAST on file(s): $INFILE "
printf "\n###################################\n"

##################
## RPSBLAST RUN ## 
##################

## FILEPATHS: I/O
rpsblast -db cdd_delta -query $INFILE -evalue 1e-5 -out ${OUTFILE} -outfmt "6 qacc sacc sseqid pident ppos length mismatch qstart qend sstart send evalue bitscore staxids"

printf "\n###################"
printf "\nEND OF RPSBLAST RUN"
printf "\n###################"
