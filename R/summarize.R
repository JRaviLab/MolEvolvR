## To summarize by lineages, DA and GC
## Modified: Jun 07, 2019
## Created: Jun 07, 2019 (was part of cleanup.R)
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzchen)

###########################
## COUNTS of DAs and GCs ##
## Before/after break up ##
###########################
## Function to obtain element counts (DA, GC)
counts_elements <- function(x, min.freq) {	x %>%
    table() %>%
    as_tibble() %>%
    `colnames<-`(c("elements", "freq")) %>%
    filter(!grepl("^-$", elements)) %>%		# remove "-"
    arrange(-freq) %>% filter(freq>=min.freq)
}

## Function to break up ELEMENTS to WORDS for DA and GC
elements2words <- function(x, type) {
  y <- x %>%
    str_replace_all("\\,"," ") %>%
    str_replace_all("\""," ")
  switch(type,
         da2doms = {z <- y %>%
           str_replace_all("\\+"," ")},
         gc2da = {z <- y %>%
           str_replace_all("\\<-"," ") %>%
           str_replace_all("-\\>"," ") %>%
           str_replace_all("\\|"," ")})
  # str_replace_all("^c\\($", " ") %>%		# remove "c("
  # str_replace_all("\\)$", " ") %>%			# remove ")"
  # str_replace_all("\\(s\\)"," ") %>%		# Ignoring repeats
  # str_replace_all("-"," ") %>%
  ## replace \n, \r, \t
  z %>%
    str_replace_all("\n"," ") %>%
    str_replace_all("\r"," ") %>%
    str_replace_all("\t"," ") %>%
    ## replace multiple spaces ...
    str_replace_all("    "," ") %>%
    str_replace_all("   "," ") %>%
    str_replace_all("  "," ") %>%
    str_replace_all("  "," ")
}

## Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
## to be used after elements2words
words2wc <- function(x){ x %>%
    str_replace_all("   "," ") %>%
    str_replace_all("  "," ") %>% str_replace_all("  "," ") %>%
    paste(collapse=" ") %>%
    strsplit(" ") %>%
    # filter(grepl(query.list[j], Query)) %>% # Create separate WCs for each Query
    # select(DA.wc) %>%
    table() %>% as_tibble() %>%
    `colnames<-`(c("words", "freq")) %>%
    ## filter out 'spurious-looking' domains
    filter(!grepl(" \\{n\\}", words)) %>%
    filter(!grepl("^c\\($", words)) %>%		# remove "c("
    filter(!grepl("^\\)$", words)) %>%		# remove ")"
    filter(!grepl("^-$", words)) %>%			# remove "-"
    filter(!grepl("^$", words)) %>%				# remove empty rows
    filter(!grepl("^\\?$", words)) %>%		# remove "?"
    filter(!grepl("^\\?\\*$", words)) %>%	# remove "?*"
    filter(!grepl("^tRNA$", words)) %>%		# remove "tRNA"
    filter(!grepl("^ncRNA$", words)) %>%	# remove "ncRNA"
    filter(!grepl("^rRNA$", words)) %>%		# remove "rRNA"
    # filter(!grepl("\\*", words)) %>%			# Remove/Keep only Query
    arrange(-freq)
}
## Function to filter based on frequencies
filter.freq <- function(x, min.freq){ x %>%
    filter(freq>=min.freq)
}

#########################
## SUMMARY FUNCTIONS ####
## Changed Lineage to Lineage.final on Aug 31
#########################
## Function to summarize and retrieve counts by Domains & Domains+Lineage
summ.DA.byLin <- function(x) { x %>%
    filter(!grepl("^-$", DomArch.norep)) %>%
    group_by(DomArch.norep, Lineage.final) %>%
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}
summ.DA <- function(x){ x %>%
    group_by(DomArch.norep) %>%
    summarise(totalcount=sum(count), totallin=n()) %>% # totallin=n_distinct(Lineage),
    arrange(desc(totallin), desc(totalcount)) %>%
    filter(!grepl(" \\{n\\}",DomArch.norep)) %>%
    filter(!grepl("^-$", DomArch.norep))
}
summ.GC.byDALin <- function(x) { x %>%
    filter(!grepl("^-$", GenContext.norep)) %>%
    filter(!grepl("^-$", DomArch.norep)) %>%
    filter(!grepl("^-$", Lineage.final)) %>% filter(!grepl("^NA$", DomArch.norep)) %>%
    group_by(GenContext.norep, DomArch.norep, Lineage.final) %>%
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}
summ.GC.byLin <- function(x) { x %>%
    filter(!grepl("^-$", GenContext.norep)) %>%
    filter(!grepl("^-$", DomArch.norep)) %>%
    filter(!grepl("^-$", Lineage.final)) %>% filter(!grepl("^NA$", DomArch.norep)) %>%
    group_by(GenContext.norep, Lineage.final) %>% # DomArch.norep,
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}
summ.GC <- function(x) { x %>%
    group_by(GenContext.norep) %>%
    summarise(totalcount=sum(count),
              totalDA=n_distinct(DomArch.norep),
              totallin=n_distinct(Lineage.final)) %>% # totallin=n_distinct(Lineage.final),
    arrange(desc(totalcount), desc(totalDA), desc(totallin)) %>%
    filter(!grepl(" \\{n\\}",GenContext.norep)) %>%
    filter(!grepl("^-$", GenContext.norep))
}


##################
#'Total Counts
#'
#'Creates a data frame with a totalcount column
#'
#'This function is designed to sum the counts column by either Genomic Context or Domain Architecture and creates a totalcount column from those sums.
#'
#' @param prot A data frame that must contain columns:
#' \itemize{\item Either 'GenContext.norep' or 'DomArch.norep' \item count}
#' @param cutoff Numeric. Cutoff for total count. Counts below cutoff value will not be shown. Default is 0.
#' @param type Character. Either "GC" for a total count by Genomic Context groupings or "DA" for a total count by Domain Architecture groupings.
#' @examples total_counts(pspa-gc_lin_counts,0,"GC")
#' @note Please refer to the source code if you have alternate file formats and/or
#' column names.
total_counts <- function(prot ,cutoff = 0, type = "GC"){
  if(type == "GC"){
    gc_count <- prot %>% group_by(GenContext.norep) %>% summarise(totalcount = sum(count))  %>% filter(totalcount >= cutoff)
    total <- left_join(prot,gc_count, by = "GenContext.norep")
  }
  else if(type == "DA"){
    da_count <- prot %>% group_by(DomArch.norep) %>% summarise(totalcount = sum(count))  %>% filter(totalcount >= cutoff)
    total <- left_join(prot,da_count, by = "DomArch.norep")
  }
  return(total)
}


#'Find Paralogs
#'
#'Creates a data frame of paralogs.
#'
#'This function returns a dataframe containing paralogs and the counts.
#'
#'@param df A data frame containing columns Species and Lineage
#'@example find_paralogs(pspa)
#'@note Please refer to the source code if you have alternate file formats and/or
#'column names.
find_paralogs <- function(df){
  #Remove eukaryotes
  df <- df %>% filter(!grepl("^eukaryota",Lineage))
  paralogTable <- select(df,Species) %>% group_by(Species)%>% count()%>% filter(n>1)
  colnames(paralogTable)[colnames(paralogTable)=="n"] = "Count"
  return(paralogTable)
}


##################################
## Descriptions for functions ####
##################################
# ## counts: Function to obtain element counts (DA, GC)
# cat("Counts DAs and GCs.
# \nFor e.g.:
# query.sub$DomArch.norep %>%
#   counts(n)
# query.sub$GenContext.norep %>%
# counts(n)")

# ## elements2words: Function to break up ELEMENTS to WORDS for DA and GC
# cat("Converting DA to domains and GC to DAs.\n2 switches: da2doms and gc2da
# \nFor e.g.:
# query.sub$DA.doms <- query.sub$DomArch.norep %>%
#   elements2words(\"da2doms\")
# query.sub$GC.da <- query.sub$GenContext.norep %>%
# 	elements2words(\"gc2da\")")


# ## words2wc: Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
# cat("Word counts for broken up domains from DAs and DAs from GCs.
# \nFor e.g.:
# DA.doms.wc <- query.sub$DA.doms %>%
#   words2wc()")