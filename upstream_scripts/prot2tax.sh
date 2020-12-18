#!/bin/bash
# Script to retrieve the TaxID for each protein

ACCNUM=$1
prot2tax()
{
	TAX=$(elink -db protein -id $ACCNUM -target taxonomy | efetch -format uid)
	printf "$ACCNUM\t$TAX\n"
}

prot2tax $ACCNUM


