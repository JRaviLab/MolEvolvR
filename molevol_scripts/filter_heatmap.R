## Currently tested w/ DcdV data

## Dependencies
library(tidyverse); library(here); library(rlang)
library(data.table)
library(d3heatmap) # https://github.com/rstudio/d3heatmap
library(heatmaply) # https://github.com/talgalili/heatmaply

source(here("scripts/colnames_molevol.R"))

## FILEPATHS
## Assuming that we are starting with "molevol.Rproj"
blast_path <- here("data/phage_defense/")
gen_path <- here("data/genomes/")
ls_path <- here("../laurensosinski/data/molevolvr_outputs/slps/")

##############
## BLAST IN ##
##############
## Starting with web BLAST outputs
## READ all BLAST files (blastx + blastp)
source_files <- list.files(paste0(blast_path, "blast", collapse=""),
                           pattern="*.txt")
# source_files <- dir(paste0(blast_path, "blast", collapse=""),
#                            pattern="*.txt")
source_files_path <- paste0(blast_path, "blast/", source_files)
web_blast_combnd <- source_files_path %>% list %>%
  pmap_dfr(fread, skip=7,
           fill=T, na.strings=c(""), header=F) %>%
  setnames(web_blastp_colnames) %>%
  filter(!is.na(Query))

## Starting with command-line BLAST outs
source_files <- dir(ls_path, pattern="^WP_.*cln.*", recursive=T)
source_files_path <- paste0(ls_path, source_files)
cl_blast_combnd <- source_files_path %>% list %>%
  pmap_dfr(fread, fill=T, na.strings=c(""), header=T)

## Not needed for cleaned up scripts
# cl_blast_combnd <- cl_blast_combnd %>%
#   transmute(Query=qacc, AccNum=sseqid,
#             Species=sscinames, TaxID=taxid, Lineage, ProtDesc=stitle,
#             PcPositive=ppos_adjusted, Evalue=evalue,
#             ClusterID, DomArch.Pfam, DomArch.Phobius) #DomArch.TMHMM

# write_tsv(cl_blast_combnd,
#           path=paste0(blast_path, "combnd.txt"),
#           col_names=T)

#################################
## Filter by RefRep, Pathogens ##
#################################
## RESTRICTING RESULTS FOR VIZ ##
## Restrict by VFDB Pathogen Genera
vfdb_genera_path <- paste0(gen_path,
                           "vfdb_pathogen_genera_list.txt")
vfdb_genera <- read_csv(vfdb_genera_path, col_names=F)

blast_sub <- cl_blast_combnd %>%
  filter(grepl(Species,
               pattern=paste0(vfdb_genera$X1, # Filter by VFDB pathogen genomes
                              collapse="|")))

## Restrict by REF+REP genomes
patric_path <- paste0(gen_path,
                      "PATRIC-firmicutes-compl_refrep_good_public_202008.csv")
patric_refrep <- read_csv(patric_path, col_names=T)
patric_colnames <- str_replace_all(string=colnames(patric_refrep),
                                   pattern=" ", replacement="_") # Rm space
colnames(patric_refrep) <- patric_colnames # Replace colnames
ref_taxIDs <- patric_refrep$NCBI_Taxon_ID  # Subset TaxIDs to filter BLAST hits

blast_sub <- cl_blast_combnd %>%
  filter(TaxID %in% ref_taxIDs) # Filter by PATRIC ref/rep genomes

##################
## DATA for VIZ ##
##################
## Pivoting long to wide
selected_col="Species"
selected_col <- sym(selected_col)

blast_sub_plot <- blast_sub         %>%
  filter(!is.na({{selected_col}}))  %>%
  group_by(Query, {{selected_col}})

blast_sub_plot_wide <- blast_sub                   %>%
  mutate(Lineage=str_replace_all(string=Lineage, pattern=">",
                                 replacement="_"),
         Species=str_replace_all(string=Species, pattern=" ",
                                 replacement="."))    %>%
  select(Query, {{selected_col}}, PcPositive)        %>%
  # filter(PcPositive>=25)                             %>% # Filter by ppos≥25
  filter(!is.na({{selected_col}}))                    %>%
  group_by(Query, {{selected_col}})                   %>%
  arrange(-PcPositive)                               %>%
  summarise(ppos_max=max(PcPositive))                %>%
  mutate(ppos_max=as.numeric(ppos_max))               %>%
  pivot_wider(names_from=Query, values_from=ppos_max) %>%
  column_to_rownames(var=as_string(selected_col))     %>%
  as.data.frame()

blast_sub_plot_wide[is.na(blast_sub_plot_wide)] = 0

blast_sub_plot_wide <- blast_sub_plot_wide %>%
  rename_at(vars(starts_with("WP_")),
            funs(str_replace(., pattern="WP_.*.1_", replacement="")))
# colnames(blast_sub_plot_wide) <- c("V.cholerae", "V.parahaemolyticus",
#                                    "A.veronii", "E.coli",
#                                    "P.mirabilis", "E.cloacae")
# "Vibrio.cholerae", "Vibrio.parahaemolyticus", "Aeromonas.veronii",
# "Escherichia.coli", "Proteus.mirabilis", "Enterobacter.cloacae"


#############
## DATAVIZ ##
#############
## HEATMAPS -- Query vs Select Spp/Lineages
## Using ggplot
ggplot(blast_sub_plot, aes(y={{selected_col}}, x=Query, fill=PcPositive)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="darkgreen") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=40, hjust = 1)) +
  labs(fill="% Similarity") +
  xlab("Query proteins") + ylab(as_string(selected_col))

## Using d3heatmap
d3heatmap(blast_sub_plot_wide, scale = "column", show_grid=F, dendrogram="row",
          colors="Blues", #k_row=8, #theme="dark",
          xaxis_font_size="8pt", yaxis_font_size="8pt",
          height="1250", width="1000")

## Using Heatmaply
heatmaply(blast_sub_plot_wide, #[,-c(1)],
          column_text_angle = 70, fontsize_row = 8, fontsize_col = 8,
          label_names=c(as_string(selected_col), "query", "ppos"),
          dendrogram="row",
          main="Homologs", key.title="% Similarity\n",
          xlab="Query proteins", #ylab="Species",
          # row_side_colors = blast_sub_plot_wide[, 1],
          # row_side_palette=ggplot2::scale_fill_gradient2(low = "white",
          #                                                high = "green",
          #                                                midpoint = 25,
          #                                                limits = c(0, 100)),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "white", high = "darkgreen",
            midpoint = 10, limits = c(0, 100)),
          branches_lwd=0.2)

