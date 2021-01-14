
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





