# CARD README

## Source:
This dataset was downloaded from the Comprehensive Antibiotic Resistance Database (CARD) in 2024-10 at https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2


CITATION:

Alcock et al. 2023. "CARD 2023: expanded curation, support for machine learning, and resistome 
prediction at the Comprehensive Antibiotic Resistance Database" Nucleic Acids Research, 
51, D690-D699. https://pubmed.ncbi.nlm.nih.gov/36263822/

## CARD SHORT NAMES

The CARD database uses standardized abbreviations, known as CARD Short Names, for AMR gene names associated with Antibiotic Resistance Ontology terms. These names are created for compatibility across data files and outputs from the Resistance Gene Identifier (RGI). Short Names for genes with 15 or fewer characters retain the original gene name, while longer names are abbreviated to uniquely represent each gene or protein. All CARD Short Names replace whitespace with underscores. For pathogen names, CARD follows the convention of capitalizing the first letter of the genus followed by the first three letters of the species in lowercase. Where applicable, CARD Short Names adopt formats such as “pathogen_gene,” “pathogen_gene_drug,” or “gene_drug.” Full lists of these abbreviations are available in the provided files:

shortname_antibiotics.tsv
shortname_pathogens.tsv"


## FASTA

The FASTA files included here contain retrieved sequences of antimicrobial resistance genes.

## Data Files Downloaded
aro_index.tsv
This file contains an index of ARO (Antibiotic Resistance Ontology) identifiers with associated GenBank accessions. Each entry includes information used to link antibiotic resistance genes to GenBank sequences.
shortname_antibiotics.tsv
Contains standardized abbreviations for antibiotics used in CARD’s short names. These abbreviations, which follow conventions from the American Society for Microbiology (ASM) and additional custom terms, provide a uniform naming system for antibiotics referenced within CARD data.

shortname_pathogens.tsv
Lists standardized abbreviations for pathogens used in CARD. Each abbreviation represents pathogen names in a condensed format, commonly the first letter of the genus followed by the first three letters of the species. This abbreviation system simplifies pathogen referencing in CARD outputs.
