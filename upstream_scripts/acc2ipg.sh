#!/bin/bash
# Script to retrieve the TaxID for each protein

ACCNUM=$1
OUTFILE=$2

acc2ipg()
{
	printf "IPG.ID\tSource\tNucAccNum\tNucStart\tNucStop\tStrand\tAccNum\tDescription\tSpecies\tSpp.Strain\tAssemblyID\n" > $OUTFILE
	epost -input $ACCNUM -db protein | efetch -format ipg | grep "RefSeq" >> $OUTFILE
}

acc2ipg $ACCNUM $OUTFILE
