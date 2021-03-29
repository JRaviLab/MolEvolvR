## Create PSP tables

## Run from molevol_scripts.Rproj
## Set relative paths accordingly
source("R/cleanup.R")
source("molevol_scripts/create_tables.R")

#####################
## PspA/Snf7 Table ##
#####################
# all <- read_csv("data/rawdata_tsv/all_clean_20200616.csv")
# query_acc <- read_tsv("data/psp/tables/acc_list_table.txt",
#                       col_names="AccNum")
all <- read_tsv("../psp_app/data/rawdata_tsv/all_clean_combined_20210329.txt")
## Doesn't exist
query_acc <- read_tsv("../psp_app/figures/tables/acc_list_table.txt",
                      col_names="AccNum") %>%
  distinct()

# query_in_all <- query_acc$AccNum[which(query_acc$AccNum %in% all$AccNum)]

## Filter by Spreadsheet list of AccNum
tables <- all %>%
  filter(AccNum %in% query_acc$AccNum) %>%
  select(AccNum, DomArch, GeneName, Species, Lineage, GenContext) %>%
  group_by(DomArch) %>%
  arrange(DomArch, GenContext)

# Subset PspA
pspa_doms <- c("PspA", "PspA(s)", "Snf7", "Snf7(s)",
               "Vps4-AAA-ATPase", "Vps4-AAA-ATPase(s)",
               "PspAA", "PspAA(s)", "PspAB", "PspAB(s)")

tables3 <- tables %>% filter_by_doms(doms_keep=pspa_doms,
                                     column="DomArch", ignore.case=T)

# Subset PspA-free
tables4 <- tables %>% filter_by_doms(doms_remove=pspa_doms,
                                     column="DomArch", ignore.case=T)

#####################
## PspA/Snf7 Table ##
#####################
tables3_gt <- create_table(df_table=tables3,
                           title="Table S3: Representative PspA/Snf7 homologs")

tables3_gt

#########################
## Non-PspA/Snf7 Table ##
#########################
tables4_gt <- create_table(df_table=tables4,
                           title="Table S4: Representative homologs of Psp cognate partner domains")
tables4_gt

# gtsave(tables4, "tables4.pdf", path=here(),
#        vwidth=400,   vheight=744,zoom =1)

#------------------------------#
##### Saving Tables as PDF ####
#------------------------------#
# webshot::install_phantomjs()

gtsave(tables3_gt, "tables3.pdf", path=here("../psp_app/figures/tables/"),
       vwidth=800, vheight=750, zoom =1)

gtsave(tables4_gt, "tables4.pdf", path=here("../psp_app/figures/tables/"),
       vwidth=800, vheight=1250, zoom =1)

library(pdftools)

pdf_convert("../psp_app/figures/tables/tables3.pdf", dpi=600)
pdf_convert("../psp_app/figures/tables/tables4.pdf", dpi=600)

## Simple table
tables3 %>%
  gt() %>%
  tab_header(
    title=md("Table S3: Representative PspA homologs"),
    subtitle="Gene, Lineage information, and Genomic Contexts grouped by Domain Architectures"
  ) %>%
  tab_source_note(md("More information available on our [webapp](https://jravilab.shinyapps.io/psp-evolution).")) %>%
  cols_align(align="right", columns=TRUE)
