# Using gt tables for Psp tables
# Created: Aug 30, 2020
# JR

library(gt)
library(here)
library(tidyverse)
library(paletteer)

source("../molevol_scripts/R/cleanup.R")
##############
## Table Fn ##
##############
## Function to create a table | theme, blues
create_table <- function(df_table=table_sub,
                         title="Table 1: Representative PspA homologs")
{
  df_table %>%
    gt() %>%
    # adding row/row groups
    # gt(rowname_col="AccNum") %>% #rowname_col = "Host", groupname_col = "PerturbagenFamily"
    tab_spanner(label = "Gene Info", columns = matches("AccNum|GeneName")) %>%
    tab_spanner(label = "Lineage Info", columns = matches("Species|Lineage")) %>%
    # tab_spanner(label = "Disease data info", columns = matches("GEO|Host|Cell|Technology")) %>%
    # with header and footer
    tab_header(
      title = md(title),
      subtitle = "Gene, Lineage info and Genomic Contexts grouped by Domain Architectures"
    ) %>%
    tab_source_note(md("More information available on our [web-app](https://jravilab.shinyapps.io/psp-evolution).")) %>%
    cols_align(align = "right", columns = TRUE) %>%
    # ## Can color data when it's numerical
    # data_color(
    #   columns = vars(Scores),
    #   colors = scales::col_numeric(
    #     palette = paletteer::paletteer_d(palette = "nord::silver_mine", direction=-1) %>%
    #       as.character(), domain = NULL),
    #   alpha = 0.8) %>%
    tab_options(
      # Headings; Titles
      heading.background.color = "#3C5488",
      heading.border.bottom.color = "#989898",
      heading.title.font.size = "12px",
      heading.subtitle.font.size = "11px",
      # Column labels
      column_labels.background.color = "#4DBBD5", #B09C85FF
      column_labels.font.size = "12px",
      # Stubs
      stub.background.color = "#4DBBD5", #B09C85FF
      stub.border.style = "dashed",
      stub.border.color = "#989898",
      stub.border.width = "1px",
      # Row groups
      row_group.background.color = "#3C5488", #FFEFDB80
      row_group.border.top.color = "#989898",
      row_group.border.bottom.style = "none",
      row_group.font.size = "12px",
      # Summary rows
      summary_row.border.color = "#989898",
      # summary_row.background.color = "#FFEBEE",
      # grand_summary_row.background.color = "#FFFFFF",
      # Table
      table.font.color = "#323232",
      table_body.hlines.color = "#989898",
      table_body.border.top.color = "#989898",
      table.font.size = "10px",
      table.width = "80%"
    )
}

#############
## Table 1 ##
#############
# all <- read_csv("data/rawdata_tsv/all_clean_20200616.csv")
# query_acc <- read_tsv("data/psp/tables/acc_list_table.txt",
#                       col_names="AccNum")
all <- read_tsv("data/rawdata_tsv/all_clean_combined_20210304.txt")
query_acc = read_tsv("figures/tables/acc_list_table.txt", col_names = "AccNum") %>%
    distinct()

# query_in_all = query_acc$AccNum[which(query_acc$AccNum %in% all$AccNum)]

## Filter by Spreadsheet list of AccNum
tables <- all %>%
  filter(AccNum %in% query_acc$AccNum) %>%
  select(AccNum, DomArch, GeneName, Species, Lineage, GenContext) %>%
  group_by(DomArch) %>%
  arrange(DomArch, GenContext)

# Subset PspA
pspa_doms = c("PspA", "PspA(s)", "Snf7", "Snf7(s)",
              "Vps4-AAA-ATPase", "Vps4-AAA-ATPase(s)",
              "PspAA", "PspAA(s)", "PspAB", "PspAB(s)")

table1 <- tables %>% filter_by_doms(doms_keep = pspa_doms,
                                    column = "DomArch", ignore.case = T)

# Subset PspA-free
table2 <- tables %>% filter_by_doms(doms_remove = pspa_doms,
                                    column = "DomArch", ignore.case = T)

#############
## Table 1 ##
#############
# Create Table 1
table1gt <- create_table(df_table=table1,
             title="Table S1: Representative PspA/Snf7 homologs")

table1gt

#############
## Table 2 ##
#############
table2gt <- create_table(df_table=table2,
             title = "Table S2: Representative PspA partner domain homologs")
table2gt

gtsave(table2gt, "table2.pdf", path = here(),
       vwidth = 400,   vheight = 744,zoom =1)

#------------------------------#
##### Saving Tables as PDF ####
#------------------------------#
# webshot::install_phantomjs()

gtsave(table1gt, "tables1.pdf", path = here("figures/tables/"),
       vwidth = 800,   vheight = 750, zoom =1)

gtsave(table2gt, "tables2.pdf", path = here("figures/tables/"),
       vwidth = 800,   vheight = 1250, zoom =1)

library(pdftools)

pdf_convert("./figures/tables/tables1.pdf", dpi = 600)
pdf_convert("./figures/tables/tables2.pdf", dpi = 600)

## Simple table
table1 %>%
  gt() %>%
  tab_header(
    title = md("Table 1: Representative PspA homologs"),
    subtitle = "Gene, Lineage info and Genomic Contexts grouped by Domain Architectures"
  ) %>%
  tab_source_note(md("More information available on our [web-app](https://jravilab.shinyapps.io/psp-evolution).")) %>%
  cols_align(align = "right", columns = TRUE)
# NRC NPG
# "#E64B35" "#4DBBD5" "#00A087" "#3C5488" "#F39B7F" "#8491B4" "#91D1C2" "#DC0000" "#7E6148" "#B09C85"