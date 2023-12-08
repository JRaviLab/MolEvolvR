#!/bin/bash

# created by Phil Calhoun on May 8th, 2019
# edited by Lo Sosinski on 2020.09.03
# standard outformat

## Load query files
QUERY_FILE=$1
PREFIX=$2
OUTDIR=$3

## iprscan able to take multiple fasta sequences in input file

printf "\n######################"
printf "\nBEGIN INTERPROSCAN RUN"
printf "\nIn files:\n$QUERY_FILE\n"
printf "\n######################\n"

## Run interproscan

## FILEPATHS: I/O
OUTFILE=$(printf "${OUTDIR}/${PREFIX}.iprscan")

## print I/O messages
printf "\nNow processing $QUERY_FILE\n"

iprscan -i ${QUERY_FILE} -b ${OUTFILE} -f TSV --cpu ${INTERPROSCAN_CPUS:-4} \
    --appl Pfam,MobiDBlite,Phobius,Coils,SignalP_GRAM_POSITIVE,SignalP_GRAM_NEGATIVE,Hamap

printf "##################\n"
printf "END OF IPRSCAN RUN\n"
printf "##################\n"
