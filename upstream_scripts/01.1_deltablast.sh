#!/bin/bash
# Script to run DELTABLAST
# Default parameters: DB=refseq; evalue=1e-5; #hits=5000
# User inputs: INFILE PREFIX OUTDIR DB NHITS EVAL
# Authors: Lauren Sosinski, Janani Ravi
# Created on: Apr 9, 2019
# Last edited: Nov 13, 2020

##################
## USAGE & HELP ##
##################

##helpFunction()
##{
##   echo ""
##   echo "Usage: $0 -i INPUT_FILE -db DB -num_hits NHITS -eval EVAL"
##   echo -e "\t-i Input file path"
##   echo -e "\t-db Database (nr or refseq; default refseq)"
##   echo -e "\t-num_hits Number of homologs included in output (5000, 10000, etc.; default 5000)"
##   echo -e "\t-eval Cut-off value (default 1e-5)"
##   echo -e "$0 -i \"/data/research/jravilab/molevol/path/to/data/fasta/*.fa\" -db \"refseq_protein\" -num_hits \"5000\" -evalue \"1e-5\""
##   exit 1 # Exit script after printing help
##}
##
##while getopts "i:db:num_hits:eval:" opt
##do
##   case "$opt" in
##      i ) INFILE="$OPTARG" ;;					# assign input file to variable
##      db ) DB="$OPTARG" ;;					# assign database to variable
##      num_hits ) N="$OPTARG" ;;				# assign number of hits
##      eval ) EVAL="$OPTARG" ;;				# assign evalue
##      h | ? ) helpFunction ;; 				# Print helpFunction in case parameter is non-existent
##      *) DB="refseq_protein", NHITS=10000, EVAL=1e-5 ;;	# default parameters for db, number of hits, and evalue
##   esac
##done

# Print helpFunction in case parameters are empty
#if [ -z "$INFILE" ]
#then
#   echo "Input file not designated";
#   helpFunction
#fi

################
## DELTABLAST ##
################

## FILEPATHS: I/O
INFILE=$1
PREFIX=$2
OUTDIR=$3
DB=$4
NHITS=$5
EVAL=$6
OUTFILE=$(printf "${OUTDIR}/${PREFIX}.dblast.tsv")

## Print I/O messages
printf "\nNow processing: $INFILE"
printf "\nRunning against: $DB\nE-value: ≤$EVAL\nTop $NHITS hits/alignments"
printf "\nOutput filepath: ${OUTFILE}\n"

## Designating different outfiles/parameters based on $DB
case "$DB" in
  nr ) dblast_db="nr" ;;
  refseq ) dblast_db="refseq_protein" ;;
esac

## core script
deltablast -query $INFILE -db $dblast_db -evalue $EVAL -num_alignments $NHITS -out "$OUTFILE" -num_threads 10 -outfmt '6 qacc sacc sseqid sallseqid stitle sscinames staxids pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore ppos'

