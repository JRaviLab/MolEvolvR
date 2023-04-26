#!/bin/bash 
# PC, LS & JR
# Last updated: Jun 9, 2020
# Created: May 8, 2019 
# To run interproscan on protein fasta files
# standard outformat

###########
## SETUP ##
###########
export INTERPRO=/mnt/research/common-data/Bio/iprscan:/mnt/research/common-data/Bio/iprscan/data:/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:/var/interproscan/data:$INTERPRO
export SIGNALP=/var/interproscan/bin/signalp/4.1

# Loading Modules and Checking Java (needs to be updated)
module purge
module use /mnt/home/johnj/software/modulefiles
module load iprscan

#########
## I/O ##
#########
#Query and Subject File
read -p "Enter INFILEPATH (e.g., /path/to/fasta/files/*.fa): " INFILE
NOW=$(date +'%Y-%m-%d %H:%M:%S')

printf "\n####################"
printf "\nRUNNING IPRSCAN"
printf "\n$NOW"
printf "\nWorking on your input files: $INFILE"
printf "\n####################\n"

###################
## INTERPROSCAN5 ## 
###################
for f in $INFILE
do
## FILEPATHS: I/O
DIR=$(dirname $f)
#FILE=$(basename -s .fa $f)
FILE=$(basename $f | sed "s/.edirect.*$//")     # when reading edirect outfile
FILE=$(basename $FILE | sed "s/.fa.*$//")       # when reading fa|faa|fasta

OUTFILE=$(printf "${DIR}/${FILE}.iprscan")

## Print I/O messages
printf "\nNow processing $f"
printf "\nOutFile: $OUTFILE"

## Core script
iprscan -i $f -b $OUTFILE -f TSV --cpu 20
done

NOW=$(date +'%Y-%m-%d %H:%M:%S')
printf "\n####################"
printf "\nEND OF IPRSCAN RUN"
printf "\n$NOW"
printf "\n####################\n"

