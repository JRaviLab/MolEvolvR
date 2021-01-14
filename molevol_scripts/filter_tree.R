
library(here)
library(tidyverse)

source(here("R/colnames_molevol.R"))
source(here("R/combine_files.R"))
source(here("R/combine_analysis.R"))
# source(here('R/ipr2viz.R'))
source(here('molevol_scripts/ipr2viz_copy.R'))


## READ FILE
inpath="../molevol_data/project_data/phage_defense/full_analysis_20210108/"

cln_combnd <- read_tsv(paste0(inpath, "cln_combined.tsv", collapse=""),
                       col_names = T) %>%
  mutate(Genus=word(Species, start=1L))

# cln <- fread(infile_full, sep ="\t", fill = T)

## Filter by DA of interest
table(cln_combnd$DomArch.Gene3D) %>%
  sort(decreasing = T)
g3d_ofinterest <- c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                    "Cytidine_Deaminase_domain_2")


## Group by Lineage, Genus, and DomArch | reverse sort by PcPositive
col1 <- "Lineage"; sym1 <- sym(col1)
col2 <- "Genus"; sym2 <- sym(col2)
col3 <- "DomArch.Pfam"; sym3 <- sym(col3)

cln <- cln_combnd %>%
  filter(grepl(g3d_ofinterest[2], DomArch.Gene3D)) %>%
  filter(grepl(g3d_ofinterest[1], DomArch.Gene3D)) %>%
  select(AccNum, Name, ProteinName,
         Query, PcPositive,
         TaxID, Species, Genus, Lineage,
         SourceDB, Completeness,
         starts_with("DomArch")) %>%
  group_by({{sym1}}, {{sym2}}, {{sym3}}) %>%
  arrange(-PcPositive) %>%
  slice_head(n=1) %>%
  filter(!is.na({{sym1}}) & !is.na({{sym2}}) & !is.na({{sym3}}))

table(cln$DomArch.Gene3D) %>%
  sort(decreasing = T)

accnum_list <- cln$AccNum %>% as.character()
write_lines(x=accnum_list,
          path=paste0(inpath, "g3d_both_accnum.txt", collapse=""))

filter_tree = function(cln_comined_path,ppos_threshold = 30,
                       subset_col = "DomArch.pfam", domains_of_interest = c(),
                       interest_col = NULL, tail_cutoff = 2){
  cln_combnd <- read_tsv(cln_comined_path,
                         col_names = T) %>%
    mutate(Genus=word(Species, start=1L)) %>%
    select(AccNum, Name, ProteinName,
           Query, PcPositive,
           TaxID, Species, Genus, Lineage,
           SourceDB, Completeness,
           starts_with("DomArch")) %>% group_by(AccNum) %>%
    arrange(-PcPositive) %>%
    slice_head(n = 1)

  # Do PPos thresholding before or after grouping?
  cln_combnd = cln_combnd %>% filter(PcPositive >= ppos_threshold)

  interest_sym =  sym(interest_col)
  for(dom in domains_of_interest)
    cln_combnd = cln_combnd %>% filter(grepl(dom, {{interest_sym}}))

  col1 <- "Lineage"; sym1 <- sym(col1)
  col2 <- "Genus"; sym2 <- sym(col2)
  col3 <- subset_col; sym3 <- sym(col3)

  hundreds = cln_combnd %>% filter(PcPositive == 100) %>% pull(AccNum)

  grp = cln_combnd %>% group_by({{sym1}},{{sym2}},{{sym3}}) %>%
    summarise(count = n(), MaxPPos = max(PcPositive)) %>% arrange(-count, -MaxPPos) %>%
    filter(!is.na({{sym3}}) & !is.na({{sym1}}))  %>% filter(count >= tail_cutoff)

  oneper = cln_combnd %>% group_by({{sym1}},{{sym2}},{{sym3}}) %>% arrange(-PcPositive) %>%
    filter(!is.na({{sym3}}) & !is.na({{sym1}})) %>%
    slice_head(n = 1)

  # For each of these, get the highest cutoff
  a = map(1:nrow(grp), function(x){
    sub = cln_combnd %>% filter(Lineage == grp$Lineage[x] &
                          Genus == grp$Genus[x] &
                            {{sym3}} == (grp %>% pull({{sym3}}))[x])
    r = sub %>% filter(PcPositive == grp$MaxPPos[x])
    return(r$AccNum[1])
  }) %>% unlist()

  res = c(a, hundreds) %>% unique()

  return(res)
}


acc = filter_tree("../molevol_data/project_data/phage_defense/full_analysis_20210108/cln_combined.tsv",
            domains_of_interest =  c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                                     "Cytidine_Deaminase_domain_2"),
            interest_col = "DomArch.Gene3D", tail_cutoff = 5)
source("R/pre-msa-tree.R")
tmp_fa = tempfile()
acc2fa(acc, tmp_fa)

tmp_msa = tempfile()

alignFasta(tmp_fa, tool = "ClustalOmega" ,outpath = tmp_msa )

source("scripts/tree.R")
tree = seq_tree(tmp_msa)


# write(read_file(tmp_msa), "../molevol_data/project_data/phage_defense/full_analysis_20210108/pfam_subset_msa.msa")
