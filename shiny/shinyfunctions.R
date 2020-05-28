library(tidyverse)
library(UpSetR)
library(gridExtra)
library(docstring)

heatmap_slider <- function(prot_lin,type){
  protcounts <- prot_lin %>% group_by(
    switch(type,
           da=DomArch.norep,
           gc=GenContext.norep
    )) %>%
    summarise(totalcount=sum(count)) %>% arrange(desc(totalcount))
  rn <- length(protcounts$totalcount)
  if(rn<=10){
    sliderinitial <- rn
  }
  else{
    sliderinitial <- 10
  }
  cutoff_init <- as.numeric(protcounts[sliderinitial,"totalcount"])
  max <- as.numeric(protcounts[1,"totalcount"])
  df <- data.frame("rn"=rn,"slider_init"=sliderinitial,
                   "cutoff_init"=cutoff_init,"max"=max)
  return(df)
}


rownumber_to_cutoff <- function(query_data, rownumber, col, RemoveAstrk = F)
{
  # Convert Rownumber input counts to percentage cutoffs that plotting functions require
  #

  # Row number of heatmap corresponds to number of GC/ DA

  if(RemoveAstrk)
  {
    query_data <- query_data %>% remove_astrk(colname = col)
  }

  column <- sym(col)

  tc <- query_data %>% total_counts(100, column = col)

  rows <- tc %>% select({{column}},'CumulativePercent') %>% distinct() #%>% arrange(CumulativePercent)

  cutoff <- rows[rownumber,'CumulativePercent'] - .001

  print(100-cutoff)

  return(100-cutoff)
}

Queries_MaxPercent <- function(prot = all, queries, colname, RemoveAstrk = F)
{
  # Summarize by the queries and lineages, and record the max percentage for each query protein
  # This allows easy conversions from Row Cutoffs to Query cutoffs
  if(RemoveAstrk)
  {
    # Remove Asterisk if necessary
    prot = remove_astrk(prot, colname)
  }

  # Total Count Summary
  prot_tc <- prot %>% total_counts(column = colname, cutoff = 100)

  colname = sym(colname)

  query_max_percs = data.frame(Query = queries, maxPercent = double(length(queries)))
  # Find largest cumulative percentage for each query
  for(x in 1:length(queries))
  {
    q = queries[x]
    q = sym(q)
    # Filter by query
    query_tc = prot_tc %>% filter(grepl(q,{{colname}}))

    # Find max Percentage of the query
    max_perc = max(query_tc$CumulativePercent)
    query_max_percs$maxPercent[x] = max_perc
  }

  query_max_percs = query_max_percs %>% arrange(-maxPercent)

  return(query_max_percs)
}



collapseInput <- function(inputId, boxId) {
  #########
  # Function to create input based off of if a box is collapsed or not

  tags$script(
    sprintf(
      "$('#%s').closest('.box').on('hidden.bs.collapse', function () {Shiny.onInputChange('%s', true);})",
      boxId, inputId
    ),
    sprintf(
      "$('#%s').closest('.box').on('shown.bs.collapse', function () {Shiny.onInputChange('%s', false);})",
      boxId, inputId
    )
  )
}
