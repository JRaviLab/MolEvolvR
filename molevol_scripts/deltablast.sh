#!/bin/bash
# JR & LS
# Created on: 2019-04-09; Edited: Apr 6, 2020

###########
## SETUP ##
###########
export BLASTDB=/mnt/research/common-data/Bio/blastdb/v5:/mnt/research/common-data/Bio/blastdb/MoreFromJanani:/mnt/research/common-data/Bio/blastdb:/mnt/research/common-data/Bio/blastdb/FASTA:$BLASTDB

## LOAD MODULES
# OPTION 1
# Loading HPC BLAST and BLAST+
#module purge
#module load BLAST/2.2.26-Linux_x86_64
#module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196 BLAST+/2.8.1-Python-2.7.14

## OPTION 2
# Loading Bioinformatics modules from JJ
# BLAST and BLAST+ are currently installed on the HPCC under 
# /mnt/home/johnj/software and /mnt/home/johnj/software/BLAST+
module use /mnt/home/johnj/software/modulefiles
module load BLAST+ BioPerl

#########
## I/O ##
#########
NOW=$(date +'%Y-%m-%d %H:%M:%S')
printf "\n##################"
printf "\nRUNNING DELTABLAST"
printf "\n$NOW"
printf "\n##################"

# Query and Subject files
#QUERY_FILE="/mnt/research/cpathogeno/molevol/data/mycobacteria/k10-fur/*.fa"
read -p "Enter INFILEPATH (e.g., /path/to/fasta/files/*.fa): " INFILE
printf "\nWorking on your input files: $INFILE \n" #DB: $DB\n"

################
## DELTABLAST ##
################
for f in $INFILE
do
## FILEPATHS: I/O
DIR=$(dirname $f)
FILE=$(basename -s .fa $f)
OUTFILE=$(printf "${DIR}/${FILE}.nr.1e5.txt")

## Print I/O messages
printf "\nNow processing: $f"
printf "\nRunning against: nr DB\nE-value: ≤1e-5\nTop 10K hits/alignments"
printf "\nOutput filepath: ${OUTFILE}\n"

## Core script
deltablast -query $f -db nr -evalue 1e-5  -num_alignments 10000 -outfmt '6 qacc sacc sseqid sallseqid stitle sscinames staxids sskingdoms pident length mismatch gapopen qstart qend sstart send evalue bitscore positive ppos slen sgi sallgi qcovs qcovhsp' -out "$OUTFILE" -num_threads 4
done

NOW=$(date +'%Y-%m-%d %H:%M:%S')
printf "\n#####################"
printf "\nEND OF DELTABLAST RUN"
printf "\n$NOW"
printf "\n#####################\n"

#########################
## UNUSED subject loop ##
#########################
# SUBJ_FILES=~/path/to/protein_fasta/ftps.protein.faa/*.faa

# for f in $SUBJ_FILES
# do 
#     echo "$f"
#     psiblast -query $QUERY_FILE -subject $f -evalue 1e-5 -num_iterations 3 -outfmt '6 sseqid sallseqid stitle staxids sblastnames sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore positive ppos slen sgi sallgi qcovs qcovhsp' -out "$f.mtb.it3.1e5.out"
# done
