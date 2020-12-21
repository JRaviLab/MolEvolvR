#!/bin/bash

# created by Lauren Sosinski 2019.01.20
# blastclust with SBATCH commands
# last edit 2020.05.18

#######################################

## Purpose/Use
# Create clusters of proteins/nucleotides of similar sequences:wq similar in function 

### NOTE: This tool is no longer maintained. In some cases, it might be more beneficial to use MMSeq2 (able to cluster up to 20,000 sequences). 

## Input: File containing FASTA sequences of proteins or nucleotides separated by identifiers
## Output: File containing one cluster per line; each sequence specified by identifier and separated by a space; sorted largest to smallest

##################
## USAGE & HELP ##
##################
# helpFunction()
# {
#   echo ""
#   echo "Usage: $0 -i INFILEPATH"
#   echo -e "\t-i input file path"
   #echo -e "\t-o Out Path"
   #echo -e "\t-db Database (NR | Refseq)"
   #echo -e "$0 -i \"/mnt/research/cpathogeno/evolvr/data/fasta/*.fsa\" -db \"nr\""
#   exit 1 # Exit script after printing help
#}

#while getopts "a:b:c:" opt
#do
#   case "$opt" in
#      i ) INFILEPATH="$OPTARG" ;;
      #db ) DB="$OPTARG" ;;
      #c ) parameterC="$OPTARG" ;;
#      h | ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
#   esac
#done

# Print helpFunction in case parameters are empty
#if [ -z "$INFILEPATH" ] #|| [ -z "$DB"]#|| [ -z "$parameterC" ]
#then
#   echo "Some or all of the parameters are empty";
#   helpFunction
#fi


################

################
## THE SCRIPT ##
################

INFILE=$1
suffix=$2
outdir=$3

printf "\n#####################################"
printf "\n## Now running BLASTCLUST on file(s): $INFILE "
printf "\n#####################################\n"

OUTFILE=$(printf "${outdir}/${suffix}.bclust.L60S80.tsv")

printf "\nPerforming BLASTCLUST analysis on $INFILE\n"
blastclust -i $INFILE -o $OUTFILE -p T -L .6 -b T -S 80 -a 8


