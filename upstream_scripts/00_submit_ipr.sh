#!/bin/bash

## COMPANION SCRIPT TO MOLEVOLVR APP ##
## Authors: Janani Ravi, Lauren Sosinski

## Input type: fasta file
## Output type: text file containing paths to individual fasta files

## USER INPUTS
INFILE=$1

dir=$(dirname $INFILE)
base=$(basename $INFILE)
cd ${dir}

find $PWD -type f -name "$base" > input.txt

INPATHS=input.txt
qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb -F $INPATHS
