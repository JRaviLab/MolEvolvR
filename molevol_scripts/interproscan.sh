#!/bin/bash 
# PC & JR
# Last updated: Apr 6, 2020; Created: May 8, 2019 
# To run interproscan on protein fasta files
# standard outformat

###########
## SETUP ##
###########
export INTERPRO=/mnt/research/common-data/Bio/iprscan:/mnt/research/common-data/Bio/iprscan/data:$INTERPRO

# Loading Modules and Checking Java (needs to be updated)
module purge
##module load Java/JDK12
##module load Java/1.12.02
##module load Java/jdk-12
module load Java/1.8.0_192
module load interproscan
module load interproscan/5.33_72.0

#########
## I/O ##
#########
#Query and Subject File
#QUERY_FILE="/mnt/research/cpathogeno/molevol/data/mycobacteria/*.fa"
read -p "Enter INFILEPATH (e.g., /path/to/fasta/files/*.fa): " INFILE
printf "\nWorking on your input files: $INFILE\n"

###################
## INTERPROSCAN5 ## 
###################
for f in $INFILE
do
## FILEPATHS: I/O
DIR=$(dirname $f)
FILE=$(basename -s .fa $f)
OUTFILE=$(printf "${DIR}/${FILE}.iprscan5")

## Print I/O messages
printf "\nNow processing $f\n"

## Core script
sh /opt/software/interproscan/5.33-72.0/interproscan.sh -i $f -b $OUTFILE
done


