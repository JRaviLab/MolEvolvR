#!/bin/bash
# Script to split a multi-fasta file to individual fasta files
# Expected user input: full FILEPATH to multi-fasta file
# Authors: Karn Jongnarangsin, Lauren Sosinski, Janani Ravi

FASTA=$1

cat $FASTA | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
        print $0 > filename
}'

#cat $FASTA | awk '{
#        if (substr($0, 1, 1)==">") {
#	  s=(substr($0,2))
#	  p=(awk '{print $1;}')
#	  filename=(p ".fa") }
#
#	print $0 > filename ;
#}'

# print filename > $(FASTA).input.txt ;
# pwd > ${FASTA}.input.txt

## look for all files in directory that end with .fa, append them to text file
# grep -wlf *.fa >> input.txt

## for loop from https://stackoverflow.com/questions/2709458/how-to-replace-spaces-in-file-names-using-a-bash-script
## renaming files -- replace all spaces in file names with underscores
# for f in *\ *; do mv "$f" "${f// /_}"; done

# awk '{print $1;}'
