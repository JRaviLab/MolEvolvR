# Script/functions to filter homologs to generate MSA and Tree

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
conflicted::conflict_prefer("filter", "dplyr")

####################################################
## Function to filter homologs to generate a tree ##
####################################################

## Example ARGS for filter_tree
inpath <- "../molevol_data/project_data/phage_defense/full_analysis_20210108/"
infile <- paste0(inpath, "cln_combined_uniq.tsv", collapse="")
cln_combined_path=paste0(inpath, "cln_combined_uniq.tsv", collapse="")
domains_of_interest=c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                      "Cytidine_Deaminase_domain_2")
subset_col1="Lineage"; subset_col2="Genus"; subset_col3="DomArch.Gene3D"
interest_col="DomArch.Gene3D"
ppos_cutoff=20 #; tail_cutoff=1

filter_tree <- function(cln_combined_path, domains_of_interest=c(),
                        subset_col1="Lineage", subset_col2="Genus",
                        subset_col3="DomArch.Pfam",
                        interest_col="DomArch.Gene3D",
                        ppos_cutoff=20, num_reps=6){ #, tail_cutoff=2
  #' Function to filter homologs to generate a tree
  #' @author Samuel Chen, Janani Ravi
  #' @param cln_combined_path Path to cln_combined.tsv (use `combine_full()`)
  #' @param domains_of_interest Default=NULL
  #' @param subset_col1 Default='Lineage'
  #' @param subset_col2 Default='Genus'
  #' @param subset_col3 Default='DomArch.Pfam'
  #' @param interest_col Default='DomArch.Gene3D'
  #' @param ppos_cutoff Default=30
  # #' @param tail_cutoff Default=2

  cln_combnd <- read_tsv(cln_combined_path, col_names=T) %>%
    mutate(Genus=word(Species, start=1L)) %>%
    select(AccNum, Name, ProteinName, Query, PcPositive,
           TaxID, Species, Genus, Lineage, SourceDB, Completeness,
           starts_with("DomArch")) %>%
    group_by(AccNum) %>%
    arrange(-PcPositive) %>%
    slice_head(n=1)

  # write_tsv(cln_combnd, col_names=T,
  #           file=paste0(inpath, "cln_combined_uniq.tsv", collapse=""))

  # Do PPos thresholding before or after grouping?
  cln_combnd <- cln_combnd %>% filter(PcPositive>=ppos_cutoff)

  interest_sym <- sym(interest_col)
  for(dom in domains_of_interest)
    cln_combnd <- cln_combnd %>% filter(grepl(dom, {{interest_sym}}))

  col1 <- subset_col1; sym1 <- sym(col1)
  col2 <- subset_col2; sym2 <- sym(col2)
  col3 <- subset_col3; sym3 <- sym(col3)

  accessions <- cln_combnd %>%
    group_by({{sym1}}, {{sym2}}, {{sym3}}) %>%
    arrange(-PcPositive) %>%
    filter(!is.na({{sym1}}) & !is.na({{sym2}}) & !is.na({{sym3}})) %>%
    slice_head(n=1) %>%
    group_by({{sym1}}, {{sym3}}) %>%
    arrange(-PcPositive) %>%
    slice_head(n=num_reps) %>% pull(AccNum)

  # groups <- cln_combnd %>%
  #   group_by({{sym1}},{{sym2}},{{sym3}}) %>%
  #   summarise(count=n(), MaxPPos=max(PcPositive)) %>%
  #   arrange(-count, -MaxPPos) %>%
  #   filter(!is.na({{sym1}}) & !is.na({{sym2}}) & !is.na({{sym3}})) %>%
  #   filter(count >= tail_cutoff)

  # oneper <- cln_combnd %>%
  #   group_by({{sym1}},{{sym2}},{{sym3}}) %>%
  #   arrange(-PcPositive) %>%
  #   filter(!is.na({{sym3}}) & !is.na({{sym1}})) %>%
  #   slice_head(n=2)

  # # For each of these, get the highest cutoff
  # accessions <- map(1:nrow(groups), function(x){
  #   sub <- cln_combnd %>%
  #     filter({{sym1}} == (groups %>% pull({{sym1}}))[x] &
  #              {{sym2}} == (groups %>% pull({{sym2}}))[x] &
  #              {{sym3}} == (groups %>% pull({{sym3}}))[x])
  #   r <- sub %>%
  #     filter(PcPositive == groups$MaxPPos[x])
  #   return(r$AccNum[1])
  # }) %>% unlist()

  # Add Query protein aka 100% identical proteins back to the AccNum list
  hundreds <- cln_combnd %>%
    filter(PcPositive == 100) %>% pull(AccNum)
  accessions <- c(accessions, hundreds) %>% unique()

  # Subset original dataframe with AccNum of interest
  cln_sub <- cln_combnd %>%
    filter(AccNum %in% accessions)

  return(accessions)
}

