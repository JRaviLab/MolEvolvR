## Functions to convert cluster files .op_ins_cls --> tsv
## Created: July 01, 2019
## Modified: Dec 11, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

## !! Note: this may not work at all since the column names here and in
## op_ins_cls, convert_opinscls_tsv.R may not match. Check files to fix issues!

#################
## Pkgs needed ##
#################
# suppressPackageStartupMessages(library(tidyverse))

##############
## COLNAMES ##
##############
## FUNCTION to ASSIGN COLUMN NAMES based on AUG 2017 VA format
# colnames.op_ins_cls <- c(
#   "AccNum", "GenContext.orig",
#   "DomArch.PFAM", "DomArch.orig", "DomArch.TMSIG",
#   "Length", "GeneName",
#   "Lineage", "Species.orig",
#   "Annotation", "GI"
# )
## FUNCTION to ASSIGN COLUMN NAMES based on DEC 2019 VA format
# colnames.op_ins_cls.clus2table <- c(
#   "AccNum", "ClustID", "ClustName.orig",
#   "GenContext.orig", "DomArch.Pfam", "DomArch.orig",
#   "-", "Length", "GeneName",
#   "Lineage", "Species.orig", "GCA_ID",
#   "Annotation", "GI"
# )

###########################
## Adding ClustIDs/Names ##
###########################
#' Clean Cluster File
#'
#' @description
#' Reads and cleans a cluster file
#'
#' This function reads a space-separated cluster file and converts it to a cleaned up data frame.
#'
#' @param path A character to the path of the cluster file to be cleaned
#' @param writepath A character designating where the tsv file of the cleaned cluster file will be written to. If value is NULL no
#' file is written. Default NULL
#' @param query A character identifying the query of the file.
#'
#' @importFrom dplyr mutate filter
#' @importFrom readr read_tsv
#' @importFrom stringr str_length
#' @importFrom tidyr separate
#'
#' @return The cleaned up cluster data frame is returned and a tsv file is written if the "writepath" parameter is used.
#'
#' @examples
#' \dontrun{
#' cleanClusterFile("data/pspa.op_ins_cls", writepath = NULL, query = "pspa")
#' }
cleanClusterFile <- function(path, writepath = NULL, query) {
    # ?? does the following line need to be changed to read_lines()?
    prot <- read_tsv(path, col_names = F)

    # ?? Unused? clust contains a column containing all the clustids
    clust <- prot %>% filter(grepl("^#", X1))

    # Separate all rows into columns by spaces and create ClustName.orig and ClustID columns
    # First warning below
    prot <- prot %>%
        separate(X1, colnames.op_ins_cls, sep = ("  +")) %>%
        mutate(ClustName.orig = "", ClustID = "")

    # ind_with_num contains a list of the row numbers with # in them.
    # This indicates that the row contains a clust id
    ind_with_num <- which(grepl("^#", prot$AccNum))

    # Separate the clustIDs (# 186;ClustName) by ";" into columns ClustID and ClustName.orig
    clsid <- separate(prot[ind_with_num, "AccNum"],
        col = AccNum,
        into = c("ClustID", "ClustName.orig"), sep = "; "
    )

    # iterate through the rows of clsid and get it into proper format
    # e.g. # 186 -> 000001.186
    for (x in 1:length(clsid$ClustID)) {
        lng <- str_length(x)
        clsid$ClustID[x] <- gsub(
            pattern = "# ",
            replacement = paste0(strrep("0", (6 - lng)), x, "."),
            clsid$ClustID[x]
        )
        # removing extra space at the start of ClustName
        clsid$ClustName.orig[x] <- gsub(
            pattern = "^ ", replacement = "",
            clsid$ClustName.orig[x]
        )
    }

    ind_with_num <- ind_with_num %>% append((length(prot$AccNum) + 1))
    # Assign CLS id and CLS name to all rows
    for (x in length(ind_with_num):2) {
        prot[which(as.numeric(rownames(prot))
        < ind_with_num[x]), "ClustName.orig"] <- clsid[x - 1, "ClustName.orig"]
        prot[which(as.numeric(rownames(prot))
        < ind_with_num[x]), "ClustID"] <- clsid[x - 1, "ClustID"]
    }

    # filter out the rows containing just the clustid and create a Query column
    prot <- prot %>%
        filter(!grepl("^#", prot$AccNum)) %>%
        mutate(Query = query)

    return(prot)
}
