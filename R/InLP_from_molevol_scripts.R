# figure 2, table 3 from the inlP paper,
# https://www.microbiologyresearch.org/docserver/fulltext/mgen/8/7/mgen000828.pdf

# fig2 code from https://github.com/JRaviLab/molevol_data/blob/main/scripts/inlp/inlp_dual_plot.R
# tbl3 code from https://github.com/JRaviLab/molevol_data/blob/main/scripts/inlp/inlp_make_table.R

library(tidyverse)
library(aplot)
library(ape)
library(ggtree)
library(tidytree)
library(seqinr)
library(msa)
library(data.table)
library(gggenes)
library(ggplot2)
library(gt)
library(here)
library(paletteer)
library(stringr)

# source("/data/research/jravilab/molevol_scripts/R/cleanup.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/summarize.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/plotting.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/networks_domarch.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/networks_gencontext.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/lineage.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/msa.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R") # in package
# source("/data/research/jravilab/molevolvr_app/scripts/ui/components.R") # doesn't exist
# source("/data/research/jravilab/molevolvr_app/scripts/MolEvolData_class.R") # no direct usage(?)
# source("/data/research/jravilab/molevolvr_app/scripts/ui/tabText.R") # no direct usage
# source("/data/research/jravilab/molevolvr_app/scripts/utility.R") # no direct usage

#' Plots the interproscan visualization aligned to a tree with inline MSA plots
#' 
#' @param app_data The app data object containing the input data.
#' @param rep The representative sequences.
#' @param phyloSelect the selected protein on the phylogeny page
#' @param tree_msa_tool The selected MSA tool on the phylogeny page (e.g. Clustal Omega, Clustal W, Muscle)
#' @param rval_phylo Indicates whether the advanced option "Phylogenetic Analysis (for homologs)" was selected for the job
#' @param acc_to_name A mapping of accessions to names.
#' @param analysis The analyses to display in the iprscan plot
#' @param show.tiplabels Whether to show tip labels on the tree plot.
#' @param progress Optional, function invoked to report progress to the caller.
#' @param verbose Optional, whether to print verbose messages, e.g. showNotification() ones
#' @return A combined plot of the phylogenetic tree and interproscan visualization.
plotDomarchTreeMSA <- function(
    app_data, rep, phyloSelect, tree_msa_tool, rval_phylo, acc_to_name,
    analysis = c("Pfam"),
    show.tiplabels = FALSE,
    progress = function(...) {}, verbose = TRUE
) {
    # ------------------------------------------------------
    # 1. read in FASTA, produce FASTA file input for MSA
    # ------------------------------------------------------
    progress(amount = 1, message = "Reading in data for FASTA file")
    
    if (!rval_phylo) {
        seqs <- readAAStringSet(app_data@fasta_path)
        names(seqs) <- sub(" .*", "", names(seqs))
        query_accession <- app_data@df %>% filter(app_data@df$QueryName == phyloSelect)
        query_accession <- unique(query_accession$Query)
        query <- seqs[query_accession]
        names(query) <- phyloSelect
        query <- AAStringSet(query)
    }

    # Generate Fasta File
    progress(amount = 1, message = "Generating FASTA file")

    rep_fasta_path <- tempfile()
    acc2fa(rep, outpath = rep_fasta_path, "sequential")
    rename_fasta(rep_fasta_path, rep_fasta_path,
        replacement_function = map_acc2name,
        acc2name = acc_to_name
    )
    if (!rval_phylo) {
        writeXStringSet(query, rep_fasta_path, append = TRUE)
    }

    # ------------------------------------------------------
    # 2. generate MSA file
    # ------------------------------------------------------
    progress(amount = 1, message = "Generating MSA file from FASTA file")

    rep_msa_path <- tempfile()
    alignFasta(rep_fasta_path, tree_msa_tool, rep_msa_path)

    # ------------------------------------------------------
    # 3. render tree from MSA file
    # ------------------------------------------------------
    progress(amount = 1, message = "Rendering tree plot from MSA")

    tree_plot <- seq_tree(
        fasta_filepath = rep_msa_path,
        show.xaxis = TRUE,
        show.tiplabels = show.tiplabels
    )

    # replace the tree plot's y axis with a new one that will
    # align properly with the iprscan plot
    tree_plot <- tree_plot +
        scale_y_continuous(
            breaks = tree_plot$data$y,
            labels = tree_plot$data$label,
            expand = c(0, 0)
        ) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank()
        )

    # ------------------------------------------------------
    # 4. pre-process iprscan data wrt. the tree
    # ------------------------------------------------------
    progress(amount = 1, message = "Preprocessing IPRScan data")

    # get all unique accessions that will end up in the iprscan plot
    unique_names <- unique(app_data@df$Name)

    # set how to filter the iprscan accessions:
    # - rep: filters to what the current query proteins are
    # - tree: filters to what is in the tree
    # - none: no filtering, just unique names
    filter_to <- "tree"

    # 'rep' contains just accessions, e.g. 'YP_009724390.1', whereas
    # the "Name" column contains a prefix presumably relating to the
    # organism, e.g. "VOrthor_Sacute_YP_009724390.1"
    # filter unique_names to entries where a corresponding entry exists in rep
    if (filter_to == "rep") {
        pattern <- paste(rep, collapse = "|")
        filtered_unique_names <- unique_names[grepl(pattern, unique_names)]
    } else if (filter_to == "tree") {
        # actually, let's get the names from tree_plot$data and then
        # show only those in the iprscan plot
        tree_labels <- trimws(tree_plot$data$label[tree_plot$data$isTip])

        filtered_unique_names <- unique_names[unique_names %in% tree_labels]
    }
    else {
        filtered_unique_names <- unique_names
    }

    # attempted replacement from domarch tab
    n <- ifelse("name" %in% colnames(data()), "name", "AccNum")

    # ------------------------------------------------------
    # 5. render the iprscan plot
    # ------------------------------------------------------

    # finally, generate the iprscan plot
    # (note that ipr2viz_web() might be needed when the job
    # is from an iprscan upload; we should unite the logic
    # from that method so we can just call ipr2viz() in all cases.)
    progress(amount = 1, message = "Rendering iprscan plot")

    interpro_plot <- ipr2viz(
        infile_ipr = app_data@ipr_path,
        infile_full = app_data@df,
        accessions = filtered_unique_names,
        tree_data = tree_plot$data,
        analysis = analysis,
        group_by = "Analysis",
        topn = NULL, # don't filter by topn, since we join to the tree anyway
        wrap.legend.text = TRUE,
        verbose = verbose, rows = 20, cols = 20
        # topn = 20, query = input$DASelect
    )


    # ------------------------------------------------------
    # 6. final combination of the two plots
    # ------------------------------------------------------

    # now let's combined them
    combined_plot <- tree_plot %>% insert_left(interpro_plot)

    return(combined_plot)
}


createDomarchTable <- function(job_dir) {
    make_table <- function() {
        df <- read_tsv(file.path(job_dir, "blast_combined.tsv"))
        df <- df %>% subset(select = c(
            QueryName, Name,
            Species, Lineage,
            AccNum, PcPositive, PcIdentity,
            DomArch.Pfam
        ))
        df <- df[order(df$QueryName), ]
        df <- df %>% mutate(DomArch.Pfam = str_replace_all(DomArch.Pfam, "\\+", "\\+<br>"))
        table <- df %>%
            gt() %>%
            fmt_markdown(columns = everything()) %>%
            cols_label(
                QueryName = "Query", Name = "Subject",
                AccNum = "Accession"
            ) %>%
            tab_source_note(md("More information available on [GitHub](https://github.com/jravilab/inlp_listeria).")) %>%
            tab_options(
                # Headings; Titles
                heading.background.color = "black",
                heading.border.bottom.color = "#989898",
                heading.title.font.size = "18px",
                heading.subtitle.font.size = "18px",
                # Column labels
                column_labels.background.color = "grey50", # B09C85FF
                column_labels.font.size = "18px",
                # Stubs
                stub.background.color = "#4DBBD5", # B09C85FF
                stub.border.style = "dashed",
                stub.border.color = "#989898",
                stub.border.width = "1px",
                # Row groups
                row_group.background.color = "#3C5488", # FFEFDB80
                row_group.border.top.color = "#989898",
                row_group.border.bottom.style = "none",
                row_group.font.size = "18px",
                # Summary rows
                summary_row.border.color = "#989898",
                # summary_row.background.color="#FFEBEE",
                # grand_summary_row.background.color="#FFFFFF",
                # Table
                table.font.color = "#323232",
                table_body.hlines.color = "#989898",
                table_body.border.top.color = "#989898",
                table.font.size = "16px",
                table.width = "90%"
            )

        return(table)

        # gtsave(table, "table.pdf")
    }

    filterByGenome <- function() {
        df <- read_tsv(file.path(job_dir, "cln_combined_no_dupes.tsv"))
        df <- df %>% add_column(Assembly = "")
        # df <- df %>% arrange(desc(PcPositive)) %>% add_column(Assembly = "") %>% group_by(Species) %>% slice(1)
        for (i in 1:nrow(df)) {
            accession <- df[i, ]$AccNum
            res <- system(paste("./find_genome.sh", accession), intern = TRUE)
            res <- str_split(res, "\t")
            tryCatch(
                {
                    assembly <- res[[1]][11]
                    print(assembly)
                    df[i, ]$Assembly <- assembly
                },
                error = function(e) {
                    print("passed")
                }
            )
        }
        df_grouped <- df %>%
            arrange(desc(PcPositive)) %>%
            group_by(Species) %>%
            slice(1)
        df <- df %>% subset(Assembly %in% df_grouped$Assembly)
        write_tsv(df, file.path(job_dir, "cln_combined_by_genome.tsv"))
    }

    filterByDomains <- function() {
        df <- read_tsv(file.path(job_dir, "cln_combined_no_dupes.tsv"))
        df <- df[order(desc(df$PcPositive)), ]
        df <- df %>%
            group_by(Species) %>%
            distinct(DomArch.Pfam, .keep_all = TRUE) %>%
            ungroup()
        write_tsv(df, file.path(job_dir, "blast_combined.tsv"))
    }

    # filterByGenome()
    # filterByDomains()
    return(make_table())
}