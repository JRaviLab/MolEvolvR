#!/bin/bash
# Script to retrieve the TaxID for each protein
# Created: Dec 16, 2020 | Janani Ravi
# Ref: https://ncbi-hackathons.github.io/EDirectCookbook/

module load edirect
# NOTE: this takes NCBI API key from the wrapper scripts

INFILE=$1	# list of accession numbers (one per line)
PREFIX=$2	# prefix to files/runs
OUTDIR=$3	# location of output files
acc2info()
{
	OUTFILE=$(printf "${OUTDIR}/${PREFIX}.acc2info.tsv")
	# Get info for the query too
	#echo ${PREFIX} >> ${INFILE}
	# print colnames
	printf "AccNum.noV\tFullAccNum\tDescription\tLength\tTaxID\tSpecies\tSourceDB\tCompleteness\n" > $OUTFILE
	# Batch input of accession numbers --> Document Summaries --> Pull necessary columns --> Output
	epost -input $INFILE -db protein | \
	efetch -format docsum | xtract -pattern DocumentSummary -def "NA" \
	-element Caption,Extra,Title,Slen,TaxId,Organism,SourceDb,Completeness >> $OUTFILE
}

acc2info_phylo()
{
	OUTFILE=$(printf "${OUTDIR}/acc2info.tsv")
	epost -input $INFILE -db protein | \
	efetch -format docsum | xtract -pattern DocumentSummary -def "NA" \
	-element Caption,Extra,Title,Slen,TaxId,Organism,SourceDb,Completeness >> $OUTFILE
}

if [ $PREFIX = "NA" ];
then
	acc2info_phylo $INFILE $PREFIX $OUTDIR
else
	acc2info $INFILE $PREFIX $OUTDIR
fi

# Unused XML alternative | parts of this could be slighly buggy
#epost -input $INFILE -db protein | \
#efetch -format xml | \
#xtract -pattern Seq-entry -block Bioseq_id -element Textseq-id_accession,Textseq-id_version -block BioSource_org -def "NA" -element Object-id_id,Org-ref_taxname,OrgName_lineage \
#-tab "\t" -sep "." -def "NA" #-element Textseq-id_accession,Textseq-id_version Object-id_id,Org-ref_taxname,OrgName_lineage
#-block BioSource_org -def "NA" -element Object-id_id,Org-ref_taxname,OrgName_lineage
