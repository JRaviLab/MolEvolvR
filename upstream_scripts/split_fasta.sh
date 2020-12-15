#!/bin/bash
# Script(s) to split multifasta file into several individual fasta files named by word #1 (aka accnum)!

INFILE=$1
dir=$(dirname $INFILE)
cd ${dir}

### ADD IF AND EVERYTHING ELSE FROM 00_split ... 

## OPTION 1
# https://unix.stackexchange.com/questions/15662/splitting-text-files-based-on-a-regular-expression
awk -F "( )|(>)" '/^>/{x=""$2".faa";}{print >x;}' $INFILE



### UNUSED ###
## OPTION 2
# https://awesomeopensource.com/project/crazyhottommy/bioinformatics-one-liners
#awk '/^>/{s=++d".fa"} {print > s}' $INFILE

## OPTION 3
# https://awesomeopensource.com/project/crazyhottommy/bioinformatics-one-liners
# https://www.golinuxcloud.com/csplit-split-command-examples-linux-unix/
# https://linuxhandbook.com/csplit-command/
#csplit -z -q -n 2 -f seq_ $INFILE /\>/ {*}  

## OPTION 4
# https://unix.stackexchange.com/questions/350523/split-files-based-on-pattern-search-split-file-name-with-pattern-we-searched
# This script will take 1 parameter as input: the target file path
#targetFile="$1"        
#targetDir=$(dirname -- "$targetFile")
#targetFile=$(basename -- "$targetFile")
#cd -P -- "$targetDir" || exit
#awk '{if(gsub(/^>/,"")){name=$0;}else{print > name".txt"}}' < "$targetFile"
