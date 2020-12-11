## Script to generate heatmaps Var1 vs Var2 for combined full analysis DFs
## Currently tested w/ DcdV and Slps data

##################
## Dependencies ##
##################
library(here)
library(tidyverse)
library(data.table)
library(rlang)
library(d3heatmap) # https://github.com/rstudio/d3heatmap
library(heatmaply) # https://github.com/talgalili/heatmaply

source(here("molevol_scripts/combine_files.R"))

################
## FILE INPUT ##
################
# inpath <- c("../molevol_data/project_data/slps/")
# cln_combnd <- read_tsv(paste0(inpath, "full_combined.tsv",
#                               collapse=""),
#                        col_names=T)

## If the combined file doesn't already exist ...
## Combining full_analysis files
inpath <- c("../molevol_data/project_data/slps/full_analysis_20201207/")
cln_combnd <- combine_files(inpath,
                            pattern="*full_analysis.txt", skip=0,
                            col_names=T)

#################################
## Filter by RefRep, Pathogens ##
#################################
gen_path <- here("../molevol_data/common_data/genomes/")

## Lineages corresp. to Gram-positive bacteria | Firmicutes + Actinobacteria
blast_sub <- cln_combnd %>%
  filter(grepl(Lineage,
               pattern=paste0(c("Actinobacteria", "Firmicutes"),
                              collapse="|")))

## PATRIC bacterial genomes | compl + ref + rep + good + public + "Firmicutes"
patric_path <- paste0(gen_path,

                      "PATRIC-firmicutes-compl_refrep_good_public_202008.csv")
## PATRIC bacterial genomes | compl + ref + rep + "human hosts"
# patric_path <- paste0(gen_path, "PATRIC_bacteria_compl_refrep_humanhost.csv")

patric_sub <- read_csv(patric_path, col_names=T)

# Cleanup colnames | Rm space
patric_colnames <- str_replace_all(string=colnames(patric_sub),
                                   pattern=" ", replacement="_")
colnames(patric_sub) <- patric_colnames

blast_sub <- blast_sub %>%
  filter(TaxID %in% patric_sub$NCBI_Taxon_ID)


## VFDB Pathogen Genera
vfdb_path <- paste0(gen_path, "vfdb_pathogen_genera_grampos.txt")
vfdb_sub <- read_csv(vfdb_path, col_names=F)

blast_sub <- blast_sub %>%
  filter(grepl(Species, pattern=paste0(vfdb_sub$X1, collapse="|")))

##################
## DATA for VIZ ##
##################
selected_col_x="Species"
selected_col_x <- sym(selected_col_x)

selected_col_y="DomArch.Pfam"
selected_col_y <- sym(selected_col_y)

# Fixing Spp, Lineage cols, subsetting by PcPositive
blast_sub_plot <- blast_sub %>%
  mutate(Lineage=str_replace_all(string=Lineage, pattern=">",
                                 replacement="_"),
         Species=str_replace_all(string=Species, pattern=" ",
                                 replacement="."))           %>%
  select({{selected_col_x}}, {{selected_col_y}}, PcPositive) %>%
  filter(PcPositive>=25)                                     %>%
  filter(grepl({{selected_col_y}}, pattern="SLH"))           %>%
  # filter(!is.na({{selected_col_y}}))  %>%
  group_by({{selected_col_x}}, {{selected_col_y}})           %>%
  arrange(-PcPositive)

## Heatmap using ggplot | Query vs Select Spp/Lineages
ggplot(blast_sub_plot, aes(y={{selected_col_x}},
                           x={{selected_col_y}},
                           fill=PcPositive)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="darkgreen") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=40, hjust = 1)) +
  labs(fill="% Similarity") +
  xlab(as_string(selected_col_x)) + ylab(as_string(selected_col_y))

## Pivoting from LONG to WIDE for other heatmap plotting functions
blast_sub_plot_wide <- blast_sub %>%
  mutate(Lineage=str_replace_all(string=Lineage, pattern=">",
                                 replacement="_"),
         Species=str_replace_all(string=Species, pattern=" ",
                                 replacement="."))    %>%
  select({{selected_col_x}}, {{selected_col_y}}, PcPositive)         %>%
  # filter(PcPositive>=50)                              %>%
  filter(!is.na({{selected_col_y}}))                  %>%
  group_by({{selected_col_x}}, {{selected_col_y}})    %>%
  arrange(-PcPositive)                                %>%
  summarise(ppos_max=max(PcPositive))                 %>%
  mutate(ppos_max=as.numeric(ppos_max))               %>%
  pivot_wider(names_from={{selected_col_x}}, values_from=ppos_max) %>%
  column_to_rownames(var=as_string(selected_col_y))   %>%
  as.data.frame()

blast_sub_plot_wide[is.na(blast_sub_plot_wide)] = 0

## Heatmap using heatmaply | Query vs Select Spp/Lineages
heatmaply(blast_sub_plot_wide, #[,-c(1)],
          column_text_angle = 70, fontsize_row = 8, fontsize_col = 8,
          label_names=c(as_string(selected_col_y),
                        as_string(selected_col_x), "ppos"),
          dendrogram="row",
          main="Homologs", key.title="% Similarity\n",
          xlab=as_string(selected_col_x), #ylab="Species",
          # row_side_colors = blast_sub_plot_wide[, 1],
          # row_side_palette=ggplot2::scale_fill_gradient2(low = "white",
          #                                                high = "green",
          #                                                midpoint = 25,
          #                                                limits = c(0, 100)),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "white", high = "darkgreen",
            midpoint = 10, limits = c(0, 100)),
          branches_lwd=0.2)

## Heatmap using d3heatmap | Query vs Select Spp/Lineages
## !! Issue: Labels getting cut !!
d3heatmap(blast_sub_plot_wide,
          scale = "column", show_grid=F, dendrogram="row",
          colors="Blues", #k_row=8, #theme="dark",
          xaxis_font_size="8pt", yaxis_font_size="8pt",
          height="1250", width="1000")
