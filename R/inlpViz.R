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

# helper function from https://github.com/JRaviLab/molevol_data/blob/main/scripts/inlp/inlp_remove_dups.R
remove_dupes <- function(job_dir) {
    df <- read_tsv(file.path(job_dir, "cln_combined_only_listeria.tsv"))
    df <- distinct(df)
    cp_df <- df
    for (row in 1:nrow(df)){
        row <- df[row,]
        dups <- subset(df, df$AccNum == row$AccNum)
        if (nrow(dups) > 1){
            dups <- dups[order(dups$PcPositive, decreasing = TRUE),]
            keeper <- dups[1,]
            dups <- dups[2:nrow(dups),]$Query
            cp_df <- subset(cp_df, !(cp_df$AccNum == keeper$AccNum & cp_df$Query %in% dups))
            print("not passed")
        }
        else{
            print("passed")
        }
    }

    write_tsv(cp_df, file.path(job_dir, "cln_combined_no_dupes.tsv"))
}


inlp_fig2 <- function(job_dir) {
    # take in a file, generate trees from that
    seq_tree <- function(fasta_filepath) {
        CPCOLS <- c(
            "#00000000", "#DDA0DD", "#EE2C2C", "#CDBE70", "#B0B099",
            "#8B2323", "#EE7600", "#EEC900", "chartreuse3", "#0000FF",
            "#FFD900", "#32CD32", "maroon4", "cornflowerblue", "darkslateblue",
            "#AB82FF", "#CD6889", "#FFA07A", "#FFFF00", "#228B22",
            "#FFFFE0", "#FFEC8B", "peru", "#668B8B", "honeydew",
            "#A020F0", "grey", "#8B4513", "#191970", "#00FF7F",
            "lemonchiffon", "#66CDAA", "#5F9EA0", "#A2CD5A", "#556B2F",
            "#EEAEEE", "thistle4", "#473C8B", "#FFB6C1", "#8B1C62",
            "#FFE4B5", "black", "#FF7F50", "#FFB90F", "#FF69B4", "#836FFF",
            "#757575", "#CD3333", "#EE7600", "#CDAD00", "#556B2F", "#7AC5CD"
        )
        my_seqs <- readAAStringSet(fasta_filepath) # , format="fasta", seek.first.rec=T)
        my_seqs_msa <- msa(my_seqs)
        my_seqs_msa_aln <- msaConvert(my_seqs_msa, type = "seqinr::alignment")

        # below was commented out, does it need to change as one of the parameters? the bottom keeps
        d <- dist.alignment(my_seqs_msa_aln, "identity")
        # as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

        ## Phylogenetic tree
        ## using package ape
        ## build neighbor-joining tree
        seqTree <- nj(d)
        # plot(seqTree, main="Phylogenetic Tree of MSA")
        groupInfo <- split(
            seqTree$tip.label,
            gsub("_\\w+", "", seqTree$tip.label)
        )
        seqTree <- groupOTU(seqTree, groupInfo)
        # ggtree(seqTree, aes(color=group),
        #        layout='circular') +
        #   geom_tiplab(size=1, aes(angle=angle))
        # https://yulab-smu.top/treedata-book/chapter4.html
        offs <- 0
        tree <- ggtree(seqTree, branch.length = "dN_vs_dS") + geom_nodelab(size = 7, col = "red") + theme_tree2()
        msaplot(tree, fasta = fasta_filepath, offset = 0.5, color = CPCOLS, bg_line = TRUE)
    }

    ipr2viz_web <- function(infile_ipr,
                            accessions,
                            analysis = c("Pfam", "Phobius", "TMHMM", "Gene3D"),
                            group_by = "Analysis", name = "Name",
                            text_size = 8, legend_name = "ShortName", cols = 5, rows = 10) {
        CPCOLS <- c(
            "#AFEEEE", "#DDA0DD", "#EE2C2C", "#CDBE70", "#B0B099",
            "#8B2323", "#EE7600", "#EEC900", "chartreuse3", "#0000FF",
            "#FFD900", "#32CD32", "maroon4", "cornflowerblue", "darkslateblue",
            "#AB82FF", "#CD6889", "#FFA07A", "#FFFF00", "#228B22",
            "#FFFFE0", "#FFEC8B", "peru", "#668B8B", "honeydew",
            "#A020F0", "grey", "#8B4513", "#191970", "#00FF7F",
            "lemonchiffon", "#66CDAA", "#5F9EA0", "#A2CD5A", "#556B2F",
            "#EEAEEE", "thistle4", "#473C8B", "#FFB6C1", "#8B1C62",
            "#FFE4B5", "black", "#FF7F50", "#FFB90F", "#FF69B4", "#836FFF",
            "#757575", "#CD3333", "#EE7600", "#CDAD00", "#556B2F", "#7AC5CD"
        )
        ## To filter by Analysis
        analysis <- paste0(analysis, collapse = "|")

        ## @SAM, colnames, merges, everything neeeds to be done now based on the
        ## combined lookup table from "common_data"
        lookup_tbl_path <- "/data/research/jravilab/common_data/cln_lookup_tbl.tsv"
        lookup_tbl <- read_tsv(lookup_tbl_path, col_names = T, col_types = lookup_table_cols)

        ## Read IPR file and subset by Accessions
        ipr_out <- read_tsv(infile_ipr, col_names = T)
        ipr_out$SignDesc <- str_replace_all(ipr_out$SignDesc, "Leucine rich repeat", "LRR")
        ipr_out$SignDesc <- str_replace_all(ipr_out$SignDesc, "Leucine Rich repeat", "LRR")
        ipr_out$SignDesc <- str_replace_all(ipr_out$SignDesc, "Leucine Rich Repeat", "LRR")
        ipr_out$ShortName <- str_replace_all(ipr_out$ShortName, "Leucine rich repeat", "LRR")
        ipr_out$ShortName <- str_replace_all(ipr_out$ShortName, "Leucine Rich repeat", "LRR")
        ipr_out$ShortName <- str_replace_all(ipr_out$ShortName, "Leucine Rich Repeat", "LRR")
        ipr_out$ShortName <- str_replace_all(ipr_out$ShortName, "LRR_4", "LRR_1")
        ipr_out$SignDesc <- str_replace_all(ipr_out$SignDesc, "LRRs \\(2 copies\\)", "LRR")
        # ipr_out <- ipr_out %>% filter(Name %in% accessions)
        ## Need to fix eventually based on 'real' gene orientation!
        ipr_out$Strand <- rep("forward", nrow(ipr_out))

        ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
        ipr_out_sub <- filter(
            ipr_out,
            grepl(pattern = analysis, x = Analysis)
        )
        # dynamic analysis labeller
        analyses <- ipr_out_sub %>%
            select(Analysis) %>%
            distinct()
        analysis_labeler <- analyses %>%
            pivot_wider(names_from = Analysis, values_from = Analysis)
        # analysis_labeler[1,] = colnames(analysis_labeler)

        ipr_out_sub$label <- paste0(" ", ipr_out_sub$Name)
        lookup_tbl <- lookup_tbl %>% select(-ShortName)
        ## @SAM, make sure the following two work with the Lookup Tables!!
        # ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by = "DB.ID")
        ## PLOTTING
        ## domains as separate arrows
        # ipr_out_sub$label <- paste0(" ", ipr_out_sub$Name)
        if (group_by == "Analysis") {
            plot <- ggplot(
                ipr_out_sub,
                aes_string(
                    xmin = 1, xmax = "SLength",
                    y = "label", label = "ShortName"
                )
            ) +
                geom_gene_arrow(fill = "white", color = "grey") +
                geom_subgene_arrow(data = ipr_out_sub, aes_string(
                    xmin = 1, xmax = "SLength", y = "label", fill = "SignDesc",
                    xsubmin = "StartLoc", xsubmax = "StopLoc"
                ), color = "NA") +
                # geom_blank(data = dummies) +
                facet_wrap(~Analysis,
                    strip.position = "top", ncol = cols,
                    labeller = as_labeller(analysis_labeler)
                ) +
                # , ncol = 1 + #scales = "free",
                scale_fill_manual(values = CPCOLS) +
                theme_minimal() +
                theme_genes2() +
                theme(
                    legend.position = "bottom",
                    legend.box = "horizontal",
                    legend.key.size = unit(0.02, "npc"),
                    legend.box.margin = margin(),
                    text = element_text(size = text_size)
                ) +
                ylab("") +
                guides(fill = guide_legend(nrow = rows))
        } else if (group_by == "Query") {
            plot <- ggplot(
                ipr_out_sub,
                aes(
                    xmin = 1, xmax = SLength,
                    y = Analysis, # y = AccNum
                    label = ShortName
                )
            ) +
                geom_subgene_arrow(data = ipr_out_sub, aes_string(
                    xmin = 1, xmax = "SLength", y = "Analysis", fill = "SignDesc",
                    xsubmin = "StartLoc", xsubmax = "StopLoc"
                ), color = "white") +
                geom_gene_arrow(fill = NA, color = "grey") +
                facet_wrap(as.formula(paste("~", name)),
                    strip.position = "top", ncol = cols,
                    labeller = as_labeller(analysis_labeler)
                ) +
                scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
                theme_minimal() +
                theme_genes2() +
                theme(
                    legend.position = "bottom",
                    legend.box = "horizontal",
                    legend.key.size = unit(0.02, "npc"),
                    legend.box.margin = margin(),
                    text = element_text(size = text_size)
                ) +
                ylab("") +
                guides(fill = guide_legend(nrow = rows))
        }
        return(plot)
    }

    theme_genes2 <- function() {
        ggplot2::theme_grey() + ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(colour = "grey80", size = 0.2),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
            axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
            strip.background = ggplot2::element_blank()
            # strip.text = ggplot2::element_blank()
        )
    }

    df <- read_tsv(file.path(job_dir, "cln_combined_by_domarch.tsv"))
    in_ipr <- read_tsv(file.path(job_dir, "ipr_combined.tsv"))

    rep_fasta_path <- tempfile()
    top_acc <- df$AccNum
    acc2fa(top_acc, outpath = rep_fasta_path, "sequential")
    rename_fasta(rep_fasta_path, rep_fasta_path,
        replacement_function = map_acc2name,
        acc2name = select(df, "AccNum", "Name")
    )

    rep_msa_path <- tempfile()
    alignFasta(rep_fasta_path, "ClustalO", rep_msa_path)

    tree <- seq_tree(rep_msa_path)
    ipr_plot <- ipr2viz_web(file.path(job_dir, "ipr_combined.tsv"), file.path(job_dir, "cln_combined_by_domarch.tsv"),
        analysis = c("Pfam"), text_size = 12
    )

    last_plot <- ipr_plot %>% insert_right(tree)

    return(last_plot)

    # ggsave("plot.png", last_plot, dpi = 600, device = "png", height = 15, width = 13)
}


inlp_tbl3 <- function(job_dir) {
    make_table <- function() {
        df <- read_tsv(file.path(job_dir, "cln_combined_by_domarch.tsv"))
        df <- df %>% subset(select = c(
            QueryName, Name,
            Species, Lineage,
            AccNum, PcPositive, PcIdentity,
            DomArch.Pfam
        ))
        df <- df[order(df$QueryName), ]
        df <- df %>% mutate(DomArch.Pfam = str_replace_all(DomArch.Pfam, "\\+", "\\+<br>"))
        df$QueryName <- str_replace_all(df$QueryName, "ivanovii", "BFirmic_Livanovii_WP_111143678.1")
        df$QueryName <- str_replace_all(df$QueryName, "monocytogenes_inLP", "BFirmic_Lmonocytogenes_WP_014601135.1")
        df$QueryName <- str_replace_all(df$QueryName, "seeligeri_1_RS12040", "BFirmic_Lseeligeri_WP_012986389.1")
        df$QueryName <- str_replace_all(df$QueryName, "seeligeri_2_RS12040", "BFirmic_Lseeligeri_WP_012986390.1")
        df$QueryName <- str_replace_all(df$QueryName, "seeligeri_3_RS12040", "BFirmic_Lseeligeri_WP_012986391.1")
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
        write_tsv(df, file.path(job_dir, "cln_combined_by_domarch.tsv"))
    }

    # filterByGenome()
    # filterByDomains()
    return(make_table())
}