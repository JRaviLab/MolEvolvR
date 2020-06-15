## To summarize by lineages, DA and GC
## Created: Jun 07, 2019
## Modified: Dec 12, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
library(rlang)
conflicted::conflict_prefer("filter", "dplyr")

###########################
## COUNTS of DAs and GCs ##
## Before/after break up ##
###########################
## Function to obtain element counts (DA, GC)
count_bycol <- function(prot=prot, column="DomArch", min.freq=1) {
  counts <- prot %>%
    select(column) %>%
    table() %>%
    as_tibble() %>%
    `colnames<-`(c(column, "freq")) %>%
    filter(!grepl("^-$", column)) %>%		# remove "-"
    filter(!is.na(column)) %>%
    arrange(-freq) %>% filter(freq>=min.freq)
  return(counts)
}

## Function to break up ELEMENTS to WORDS for DA and GC
elements2words <- function(prot, column= "DomArch",conversion_type="da2doms") {
  z1 <- prot %>%
    select(column) %>%
    str_replace_all("\\,"," ") %>%
    str_replace_all("\""," ")
  switch(conversion_type,
         da2doms = { z2 <- z1 %>%
           str_replace_all("\\+"," ")},
         gc2da = { z2 <- z1 %>%
           str_replace_all("\\<-"," ") %>%
           str_replace_all("-\\>"," ") %>%
           str_replace_all("\\|"," ")})
  # str_replace_all("^c\\($", " ") %>%		# remove "c("
  # str_replace_all("\\)$", " ") %>%			# remove ")"
  # str_replace_all("\\(s\\)"," ") %>%		# Ignoring repeats
  # str_replace_all("-"," ") %>%
  ## replace \n, \r, \t
  z3 <- z2 %>%
    str_replace_all("\n"," ") %>%
    str_replace_all("\r"," ") %>%
    str_replace_all("\t"," ") %>%
    ## replace multiple spaces ...
    str_replace_all("    "," ") %>%
    str_replace_all("   "," ") %>%
    str_replace_all("  "," ") %>%
    str_replace_all("  "," ")
  return(z3)
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
    filter(!grepl("^X$|^X\\(s\\)$", words)) %>% # remove "X" and "X(s)"

    # filter(!grepl("\\*", words)) %>%			# Remove/Keep only Query
    arrange(-freq)
}
## Function to filter based on frequencies
filter_freq <- function(x, min.freq){ x %>%
    filter(freq>=min.freq)
}

#########################
## SUMMARY FUNCTIONS ####
#########################
summarize_bylin <- function(prot="prot", column="DomArch", by="Lineage",
                            query="PspA"){
  column <- sym(column); by <- sym(by)
  if(query== "all"){
    prot <- prot
  } else { prot <- prot %>% filter(grepl(pattern=query, x={{column}},
                                         ignore.case=T))}
  prot %>%
    filter(!grepl("^-$", {{column}})) %>%
    group_by({{column}}, {{by}}) %>%
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}


## Function to summarize and retrieve counts by Domains & Domains+Lineage
summ.DA.byLin <- function(x) { x %>%
    filter(!grepl("^-$", DomArch)) %>%
    group_by(DomArch, Lineage) %>%
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}

## Function to retrieve counts of how many lineages a DomArch appears in
summ.DA <- function(x){ x %>%
    group_by(DomArch) %>%
    summarise(totalcount=sum(count), totallin=n()) %>% # totallin=n_distinct(Lineage),
    arrange(desc(totallin), desc(totalcount)) %>%
    filter(!grepl(" \\{n\\}",DomArch)) %>%
    filter(!grepl("^-$", DomArch))
}
summ.GC.byDALin <- function(x) { x %>%
    filter(!grepl("^-$", GenContext)) %>%
    filter(!grepl("^-$", DomArch)) %>%
    filter(!grepl("^-$", Lineage)) %>% filter(!grepl("^NA$", DomArch)) %>%
    group_by(GenContext, DomArch, Lineage) %>%
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}
summ.GC.byLin <- function(x) { x %>%
    filter(!grepl("^-$", GenContext)) %>%
    filter(!grepl("^-$", DomArch)) %>%
    filter(!grepl("^-$", Lineage)) %>% filter(!grepl("^NA$", DomArch)) %>%
    group_by(GenContext, Lineage) %>% # DomArch.norep,
    summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
    arrange(desc(count))
}
summ.GC <- function(x) { x %>%
    group_by(GenContext) %>%
    summarise(totalcount=sum(count),
              totalDA=n_distinct(DomArch),
              totallin=n_distinct(Lineage)) %>% # totallin=n_distinct(Lineage),
    arrange(desc(totalcount), desc(totalDA), desc(totallin)) %>%
    filter(!grepl(" \\{n\\}",GenContext)) %>%
    filter(!grepl("^-$", GenContext))
}


##################
total_counts <- function(prot, column = "DomArch",
                         cutoff = 90, RowsCutoff = FALSE, digits = 2
                         #type = "GC"
){
  #'Total Counts
  #'
  #'Creates a data frame with a totalcount column
  #'
  #'This function is designed to sum the counts column by either Genomic Context or Domain Architecture and creates a totalcount column from those sums.
  #'
  #' @param prot A data frame that must contain columns:
  #' \itemize{\item Either 'GenContext' or 'DomArch.norep' \item count}
  #' @param cutoff Numeric. Cutoff for total count. Counts below cutoff value will not be shown. Default is 0.
  #' @param column Character. The column to summarize
  #' @examples total_counts(pspa-gc_lin_counts,0,"GC")
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  prot <- summarize_bylin(prot, column, by = "Lineage", query = "all")
  column <- sym(column)
  col_count <-  prot %>% group_by({{column}}) %>% summarise(totalcount = sum(count))

  total <- left_join(prot,col_count, by = as_string(column))

  sum_count <- sum(total$count)
  total <- total %>% mutate("IndividualCountPercent" = totalcount/sum_count*100) %>%  arrange(-totalcount,-count)

  cumm_percent <- total %>% select({{column}}, totalcount) %>% distinct() %>% mutate("CumulativePercent"=0)
  total_counter = 0
  for(x in length(cumm_percent$totalcount):1){
    total_counter = total_counter + cumm_percent$totalcount[x]
    cumm_percent$CumulativePercent[x] = total_counter/sum_count * 100
  }

  cumm_percent <- cumm_percent %>% select( CumulativePercent, {{column}})

  total <- total %>% left_join(cumm_percent, by = as_string(column))

  # Round the percentage columns
  total$CumulativePercent <- total$CumulativePercent %>% round(digits = digits)
  total$IndividualCountPercent <- total$IndividualCountPercent %>% round(digits = digits)

  if(RowsCutoff)
  {
    # If total counts is being used for plotting based on number of rows,
    # don't include other observations that fall below the cummulative percent cutoff
    #, but that have the same 'totalcount' number as the cutoff observation
    total <- total %>% filter(CumulativePercent >= 100-cutoff-.0001)
    return(total)
  }

  # Include observations that fall below the cummulative percent cutoff,
  # but that have the same 'totalcount' as the cutoff observation
  t <- total %>% filter(CumulativePercent >= 100-cutoff)
  if(length(t) == 0){
    cutoff_count = 0
  }
  else{
    cutoff_count = t$totalcount[nrow(t)]
  }

  total <- total %>% filter(totalcount >= cutoff_count) %>% ungroup()

  return(total)
}



# total_counts_by_query <- function(query_data, queries, colname,cutoff, RemoveAstrk = F)
# {
#   ## Get the total counts by the Queries.
#
#   lineage_by_query <- function(data, query, column, by){
#     # Function to filter data by a query and summarise it by lineage
#
#     column <- sym(column); by <- sym(by)
#
#     # filter the protein by the query
#     data <- data %>% filter(grepl(pattern=query, x={{column}},
#                                   ignore.case=T)) %>% select({{by}})
#
#     data$Query <- query
#
#     data <- data %>% filter(!grepl("^-$", {{by}})) %>%
#       group_by(Query, {{by}}) %>%
#       summarise(count=n()) %>%
#       arrange(desc(count))
#     return(data)
#   }
#
#
#   col <- sym(colname)
#
#   # query_data contains all rows that possess a lineage
#   query_data <- query_data %>% filter(grepl("a", Lineage))
#
#   query_data <- shorten_lineage(query_data, "Lineage")
#
#   if(RemoveAstrk){
#     # Remove Asterisk (*) if necessary
#     query_data <- query_data %>% remove_astrk(colname = colname)
#   }
#
#   query_lin_counts = data.frame("Query" = character(0), "Lineage" = character(0), "count"= integer())
#
#   for(q in queries){
#     # iterate over queries and do filtering/Lineage by query
#     query_lin <- lineage_by_query(data = query_data, query = q, column = {{col}}, by = "Lineage")
#     query_lin_counts <- dplyr::union(query_lin_counts, query_lin)
#   }
#
#   # Total counts of each lineage
#   lin_count_totals <- query_lin_counts %>%group_by(Lineage) %>% summarize(total_c = sum(count)) %>% arrange(total_c)
#   sum_lin <- sum(lin_count_totals$total_c)
#
#   # Calculate the percent each lineage makes up
#   lin_count_totals$IndividualPercent = lin_count_totals$total_c/sum_lin * 100
#   lin_count_totals <- lin_count_totals %>% mutate("CumulativePercent"=0)
#   total_counter = 0
#   for(x in 1:length(lin_count_totals$IndividualPercent)){
#     total_counter = total_counter + lin_count_totals$IndividualPercent[x]
#     lin_count_totals$CumulativePercent[x] = total_counter
#   }
#
#   query_total_counts <- left_join(query_lin_counts, lin_count_totals, by = "Lineage")
#
#   query_total_counts <- query_total_counts%>% group_by(total_c,count) %>% arrange(-total_c, -count)
#   # Get lineages that are above the cutoff percentage value
#   count_cutoff <- (query_total_counts %>% filter(CumulativePercent >= (100-cutoff)))
#
#   query_total_counts <- (query_total_counts %>% filter(total_c >= count_cutoff$total_c[nrow(count_cutoff)]))
#
#   # Round to 3 digits
#   query_total_counts$CumulativePercent <- query_total_counts$CumulativePercent %>% round(digits = 3)
#   query_total_counts$IndividualPercent <- query_total_counts$IndividualPercent %>% round(digits = 3)
#
#   return(query_total_counts)
# }






find_paralogs <- function(prot){
  #'Find Paralogs
  #'
  #'Creates a data frame of paralogs.
  #'
  #'This function returns a dataframe containing paralogs and the counts.
  #'
  #'@param prot A data frame filtered by a Query, containing columns Species and Lineage
  #'@example find_paralogs(pspa)
  #'@note Please refer to the source code if you have alternate file formats and/or
  #'column names.


  #Remove eukaryotes
  prot <- prot %>% filter(!grepl("^eukaryota",Lineage))
  paralogTable <- prot %>% group_by(GCA_ID) %>%
    count(DomArch)%>%
    filter(n>1 & GCA_ID != "-") %>% arrange(-n) %>% distinct()
    #%>% ccNums" = filter(prot, grepl(GCA_ID, prot))$AccNum)
  paralogTable$AccNums <- map(paralogTable$GCA_ID, function(x) filter(prot,grepl(x,GCA_ID))$AccNum)
  colnames(paralogTable)[colnames(paralogTable)=="n"] = "Count"
  ###Merge with columns: AccNum,TaxID, and GCA/ Species?
  paralogTable <- prot %>% select(Species , GCA_ID, Lineage) %>%
    right_join(paralogTable, by= c("GCA_ID")) %>%
    filter(!is.na(Count)) %>% arrange(-Count) %>% select(-GCA_ID) %>% distinct()
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
# query.sub$GenContext %>%
# counts(n)")

# ## elements2words: Function to break up ELEMENTS to WORDS for DA and GC
# cat("Converting DA to domains and GC to DAs.\n2 switches: da2doms and gc2da
# \nFor e.g.:
# query.sub$DA.doms <- query.sub$DomArch.norep %>%
#   elements2words(\"da2doms\")
# query.sub$GC.da <- query.sub$GenContext %>%
# 	elements2words(\"gc2da\")")


# ## words2wc: Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
# cat("Word counts for broken up domains from DAs and DAs from GCs.
# \nFor e.g.:
# DA.doms.wc <- query.sub$DA.doms %>%
#   words2wc()")