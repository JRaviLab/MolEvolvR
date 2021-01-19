########################################################################
## Example function calls in action: Full combined -> Fasta + MSA + Tree
## cln_full -> filtered AccNum -> fa -> cln fa -> msa -> tree -> pdf  ##
########################################################################

source("./R/pre-msa-tree.R")
source("./molevol_scripts/filter_tree.R")
source("scripts/tree.R")

## Filter homologs by Analysis | Domains of interest | PcPositive | Tail cutoffs
inpath <- "../molevol_data/project_data/phage_defense/full_analysis_20210108/"
infile <- paste0(inpath, "cln_combined.tsv", collapse="")
accessions <- filter_tree(cln_combined_path=infile,
                          domains_of_interest=c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                                                "Cytidine_Deaminase_domain_2"),
                          subset_col1="Lineage", subset_col2="Genus",
                          subset_col3="DomArch.Gene3D",
                          interest_col="DomArch.Gene3D",
                          ppos_cutoff=20, tail_cutoff=1)

## Generate Fasta file
tmp_fa <- tempfile()
acc2fa(accessions, tmp_fa)

## Cleanup Fasta file with 'Name' instead of Accessions + Annotations
cln_combined <- fread(infile, sep="\t", fill=T)

renamed <- rename_fasta(fa_path=tmp_fa,
                        outpath=tmp_fa, replacement_function=map_acc2name,
                        acc2name=cln_combined)
view(renamed)

# out_fa <- gsub(x=infile, pattern="cln_combined.tsv", replacement="cln_sub.fa")
# write(x=read_file(tmp_fa), file=out_fa)

## Generate MSA
tmp_msa <- tempfile()
alignFasta(tmp_fa, tool="ClustalO", outpath=tmp_msa)

# out_msa <- gsub(x=infile, pattern="cln_combined.tsv",
#                 replacement="cln_sub_clustalo.msa")
# write(x=read_file(tmp_msa), file=out_msa)


## Generate Tree
tree <- seq_tree(tmp_msa)
tree

