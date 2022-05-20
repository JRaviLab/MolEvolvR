#!/bin/bash 
# Script to convert AccNum to Fasta | using edirect
# Created: ??
# Authors: Lauren Sosinski, Janani Ravi
# Last modified: Aug 22, 2020

# REF: NCBI e-utilities book https://www.ncbi.nlm.nih.gov/books/NBK179288/


printf "\n####################\n"
printf "BEGIN EDIRECT SEARCH\n"
printf "####################\n"

INFILE=$1
PREFIX=$2
OUTDIR=$3
OUTFILE=$(printf "${OUTDIR}/${PREFIX}.all_accnums.fa")		# creating the output file name

printf "$INFILE\n"

## print statement for current step
printf "\nCreating temporary files\n"
cols=$(head -n 1 $INFILE | awk '{print NF}')
if [ $cols -gt 1 ]
then
# taking 2nd column in blast input file and taking only 1 copy of each homolog accession,
cat ${INFILE} | awk -F "\t" '{ print $2 }' | sort -u >> ${OUTDIR}/${PREFIX}.all_accnums.txt
else
cat ${INFILE} | awk -F "\t" '{ print $1 }' | sort -u >> ${OUTDIR}/${PREFIX}.all_accnums.txt
fi
split -l 1000 -e ${OUTDIR}/${PREFIX}.all_accnums.txt ${OUTDIR}/acc 	# split accessions up into files of 1000 accession numbers

## print statement for current step
printf "\nObtaining fasta files\n"
for x in ${OUTDIR}/acc*						# looping through each file following the pattern
do
   #sleep 1
   accnum=$(cat $x | tr '\n' ',')				# separate each accession by comma instead of newline, required for edirect search
   efetch -db protein -format fasta -id "$accnum" >> $OUTFILE	# fetch each accession's fasta file from NCBI database
done
num_seqs=$(grep ">" ${OUTFILE} | wc -l | grep -Eo "^[[:digit:]]+")
if [ "$num_seqs" -lt 1 ]
then
for x in ${OUTDIR}/acc*
do
accnum=$(cat $x | tr '\n' ',')
curl -X GET --header 'Accept:text/x-fasta' 'https://www.ebi.ac.uk/proteins/api/proteins?accession='"${accnum}" >>$OUTFILE 
done
fi
## print statement for current step
printf "\nRemoving temporary files\n"
rm ${OUTDIR}/acc*							# removes each file created by "split" function above

printf "#####################\n"
printf "END OF EDIRECT SEARCH\n"
printf "#####################\n"

