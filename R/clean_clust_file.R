library(tidyverse)

##############
## COLNAMES ##
##############
## FUNCTION to ASSIGN COLUMN NAMES based on AUG 2017 VA format
colnames.op_ins_cls <- c("AccNum", "GenContext.orig",
                         "arch.PFAM", "DomArch.orig", "arch.TMSIG",
                         "Length", "GenName",
                         "Lineage", "Species.orig",
                         "Annotation", "GI")

###########################
## Adding ClustIDs/Names ##
###########################
#'Clean Cluster File
#'
#'Reads and cleans a cluster file
#'
#'This function reads a space-separated cluster file and converts it to a cleaned up data frame.
#'The cleaned up cluster data frame is returned and a tsv file is written if the "writepath" parameter is used.
#'
#' @param path A character to the path of the cluster file to be cleaned
#' @param writepath A character designating where the tsv file of the cleaned cluster file will be written to. If value is NULL no
#' file is written. Default NULL
#' @param query A character identifying the query of the file.
#' @examples clean_clust_file("data/pspa.op_ins_cls", writepath=NULL, query="pspa")

clean_clust_file <- function(path, writepath=NULL, query){
  # ?? does the following line need to be changed to read_lines()?
  prot <- read_tsv(path, col_names=F)

  # ?? Unused? clust contains a column containing all the clustids
  clust <- prot %>% filter(grepl("^#", X1))

  #Separate all rows into columns by spaces and create ClustName and ClustID columns
  #First warning below
  prot <- prot %>%
    separate(X1, colnames.op_ins_cls, sep=("  +")) %>%
    mutate(ClustName="", ClustID="")

  #ind_with_num contains a list of the row numbers with # in them.
  #This indicates that the row contains a clust id
  ind_with_num <- which(grepl("^#", prot$AccNum))

  #Separate the clustIDs (# 186;ClustName) by ";" into columns ClustID and ClusName
  clsid <- separate(prot[ind_with_num, "AccNum"], col=AccNum,
                    into=c("ClustID", "ClusName"), sep="; ")

  #iterate through the rows of clsid and get it into proper format
  #e.g. # 186 -> 000001.186
  for(x in 1:length(clsid$ClustID)){
    lng <- str_length(x)
    clsid$ClustID[x] <- gsub(pattern="# ",
                             replacement=paste0(strrep("0", (6-lng)), x, "."),
                             clsid$ClustID[x])
    # removing extra space at the start of ClustName
    clsid$ClusName[x] <- gsub(pattern="^ ", replacement="", clsid$ClusName[x])
  }

  ind_with_num <- ind_with_num %>% append((length(prot$AccNum)+1))
  #Assign CLS id and CLS name to all rows
  for (x in length(ind_with_num):2) {
    prot[which(as.numeric(rownames(prot))
               < ind_with_num[x]), "ClustName"]= clsid[x-1, "ClusName"]
    prot[which(as.numeric(rownames(prot))
               < ind_with_num[x]), "ClustID"]= clsid[x-1, "ClustID"]
  }

  #filter out the rows containing just the clustid and create a Query column
  prot <- prot %>%
    filter(!grepl("^#", prot$AccNum)) %>%
    mutate(Query=query)

  return(prot)
}