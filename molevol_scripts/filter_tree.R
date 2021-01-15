
library(here)
library(tidyverse)

####################################################
## Function to filter homologs to generate a tree ##
####################################################
filter_tree <- function(cln_combined_path, domains_of_interest=c(),
                        subset_col="DomArch.Gene3D",
                        interest_col="DomArch.Gene3D",
                        ppos_cutoff=20, tail_cutoff=2){
  #' Function to filter homologs to generate a tree
  #' @author Samuel Chen, Janani Ravi
  #' @param cln_combined_path Path to cln_combined.tsv (use `combine_full()`)
  #' @param domains_of_interest Default=NULL
  #' @param subset_col Default='DomArch.Gene3D'
  #' @param interest_col Default='DomArch.Gene3D'
  #' @param ppos_cutoff Default=30
  #' @param tail_cutoff Default=2

  cln_combnd <- read_tsv(cln_combined_path, col_names=T) %>%
    mutate(Genus=word(Species, start=1L)) %>%
    select(AccNum, Name, ProteinName, Query, PcPositive,
           TaxID, Species, Genus, Lineage, SourceDB, Completeness,
           starts_with("DomArch")) %>%
    group_by(AccNum) %>%
    arrange(-PcPositive) %>%
    slice_head(n=1)

  # Do PPos thresholding before or after grouping?
  cln_combnd <- cln_combnd %>% filter(PcPositive>=ppos_cutoff)

  interest_sym <- sym(interest_col)
  for(dom in domains_of_interest)
    cln_combnd <- cln_combnd %>% filter(grepl(dom, {{interest_sym}}))

  col1 <- "Lineage"; sym1 <- sym(col1)
  col2 <- "Genus"; sym2 <- sym(col2)
  col3 <- subset_col; sym3 <- sym(col3)

  hundreds <- cln_combnd %>%
    filter(PcPositive == 100) %>% pull(AccNum)

  groups <- cln_combnd %>%
    group_by({{sym1}},{{sym2}},{{sym3}}) %>%
    summarise(count=n(), MaxPPos=max(PcPositive)) %>%
    arrange(-count, -MaxPPos) %>%
    filter(!is.na({{sym3}}) & !is.na({{sym1}})) %>%
    filter(count >= tail_cutoff)

  oneper <- cln_combnd %>%
    group_by({{sym1}},{{sym2}},{{sym3}}) %>%
    arrange(-PcPositive) %>%
    filter(!is.na({{sym3}}) & !is.na({{sym1}})) %>%
    slice_head(n=1)

  # For each of these, get the highest cutoff
  accessions <- map(1:nrow(groups), function(x){
    sub <- cln_combnd %>% filter(Lineage == groups$Lineage[x] &
                                   Genus == groups$Genus[x] &
                                   {{sym3}} == (groups %>% pull({{sym3}}))[x])
    r <- sub %>% filter(PcPositive == groups$MaxPPos[x])
    return(r$AccNum[1])
  }) %>% unlist()

  # Add Query protein aka 100% identical proteins back to the AccNum list
  accessions <- c(accessions, hundreds) %>% unique()

  # Subset original dataframe with AccNum of interest
  cln_sub <- cln_combnd %>%
    filter(AccNum %in% accessions)

  return(accessions)
}

## Filter homologs
inpath <- "../molevol_data/project_data/phage_defense/full_analysis_20210108/"
acc <- filter_tree(cln_combined_path=paste0(inpath, "cln_combined.tsv", collapse=""),
                   domains_of_interest=c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                                         "Cytidine_Deaminase_domain_2"),
                   subset_col="DomArch.Gene3D", interest_col="DomArch.Gene3D",
                   ppos_cutoff=20, tail_cutoff=2)

cln_sub <- filter_tree(cln_combined_path=paste0(inpath, "cln_combined.tsv", collapse=""),
                       domains_of_interest=c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                                             "Cytidine_Deaminase_domain_2"),
                       subset_col="DomArch.Gene3D", interest_col="DomArch.Gene3D",
                       ppos_cutoff=20, tail_cutoff=2)

## Generate MSA
source("R/pre-msa-tree.R")
tmp_fa <- tempfile()
acc2fa(accessions, tmp_fa)

tmp_msa <- tempfile()

alignFasta(tmp_fa, tool="ClustalO", outpath=tmp_msa)

## Generate Tree
source("scripts/tree.R")
tree <- seq_tree(tmp_msa)
tree

# write(read_file(tmp_msa), "../molevol_data/project_data/phage_defense/full_analysis_20210108/pfam_subset_msa.msa")
