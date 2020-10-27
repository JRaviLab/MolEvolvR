### Function to find the top n rows/groupings: Used to find cutoff Vals
library(tidyverse)



# top_n_cutoff <- function(prot, column, n_rows, RemoveAstrk = FALSE){
#   # column <- sym(column)
#   if(RemoveAstrk)
#   {
#     prot <- prot %>% remove_astrk(colname = column)
#   }
#
#   # Summarize by lineage and get total counts
#   prot_count <- total_counts(prot, column, cutoff = 100) %>%
#     filter(count != 1) # Remove singles
#
#   row_len <- nrow(prot_count)   # Number of rows corresponds to the amount of tiles
#   # Restrict by columns of heatmap and rows of heatmap separately?
#
#   if(row_len <= n_rows){
#     # Number of rows is less than the max rows
#     cutoff_perc <- floor(prot_count$CummulativePercent[row_len])
#   }
#   else{
#     # Number of rows exceeds limit and cutoff must be established
#     cutoff_perc <- floor(prot_count$CummulativePercent[n_rows])
#   }
#
#   return(cutoff_perc)
# }

top_n_rows_cutoff <- function(prot, column, n_rows){
  # Obtain counts of all DomArchs/GC/Lineages
  column <- sym(column)
  prot_counts <- prot %>%
    group_by({{column}}) %>%
    summarise(count = n()) %>%
    #filter(count != 1) %>%  # Filter out singles
    arrange(desc(count))

  # Sum of counts
  total_count <- sum(prot_counts$count)

  # Individual percentages
  prot_counts <- prot_counts %>%
    mutate("IndPercentage" = count / total_count * 100) %>%
    mutate("CummPercentage" = 0)

  # Rows is the total number of unique GCs/DAs/Lineages
  row_len <- nrow(prot_counts)

  # print(paste0("Row Length: ", row_len))
  # Cumulative percentages
  total_counter = 0
  for(i in row_len:1){
    total_counter = total_counter + prot_counts$count[i]
    prot_counts$CummPercentage[i] = total_counter/total_count * 100
  }

  prot_counts <- prot_counts %>% filter(!grepl("^-$", {{column}}))

  if(row_len < n_rows){
    # Less than limit, can include all data
    cutoff_perc <- 100
  }
  else{
    # Exceeds row limit
    cutoff_perc <- 100 -prot_counts$CummPercentage[n_rows]
  }
  return(cutoff_perc)

}

RowNums <- function(prot, column){
  column <- sym(column)
  prot_counts <- prot %>%
    group_by({{column}}) %>%
    summarise(count = n()) %>%
    #filter(count != 1) %>%  # Filter out singles
    arrange(desc(count))

  # Rows is the total number of unique GCs/DAs/Lineages
  row_len <- nrow(prot_counts)
  return(row_len)
}



max_word_percents <- function(data, colname)
{
  # Get The max rows (words) of the upset plots
  # Get the max percent for each word

  # Get Total Counts
  # colname = string(colname)
  tc <- data %>% total_counts(column =  colname, cutoff = 100, RowsCutoff = F, digits = 5)

  column <- sym(colname)

  if(grepl("GenContext", colname, ignore.case = T))
  {
    type = "gc2da"
  }
  else if(grepl("DomArch|ClustName", colname, ignore.case = T))
  {
    type = "da2doms"
  }

  elements_Percs = tc %>% select({{column}}, CumulativePercent) %>% distinct()
  # Get words from filter
  words.tc <- elements_Percs %>% select({{column}}) %>% elements2words(column = colname,conversion_type = type)
  # names(words.tc)[1] <- "words"
  words.tc <- words.tc %>% str_split(pattern = " ") %>% unlist()
  words.tc <- words.tc[2:(length(words.tc)-1)]

  max_perc_words = data.frame("Words" = words.tc, "MaxPercent" = double(length(words.tc)), stringsAsFactors = F)
  max_perc_words = max_perc_words %>% distinct()

  listwords = max_perc_words$Words

  listwords = str_replace_all(listwords, "\\(", "\\\\(") %>%
    str_replace_all("\\)", "\\\\)") %>%
    str_replace_all("\\+", "\\\\+")
  for(i in 1:nrow(max_perc_words))
  {
    # Iterate over all words and find their max percent
    cur_word = listwords[i]
    # Filter by the word
    word_percs = elements_Percs %>% filter(grepl(pattern = cur_word,x = {{column}}))

    max_perc = max(word_percs$CumulativePercent)

    max_perc_words$MaxPercent[i] = max_perc
  }

  max_perc_words = max_perc_words %>% arrange(-MaxPercent)
  return(max_perc_words)
}



# GetLineageCounts <- function(prot,lineage_col = "Lineage")
# {
#   # Get the row numbers for all by using the total counts
#   # Also get the initial count percentage,
#   lineage_col <- sym(lineage_col)
#
#   all_lins <- prot %>% group_by({{lineage_col}}) %>% summarise(count = n()) %>% arrange(-count)
#
#
#   sum_count <- sum(all_lins$count)
#   all_lins <- all_lins %>% mutate("IndividualPercent" = count/sum_count*100) %>%
#     mutate("CumulativePercent" = 0)
#
#   total_counter = 0
#   for(x in nrow(all_lins):1){
#     total_counter = total_counter + all_lins$count[x]
#     all_lins$CumulativePercent[x] = total_counter/sum_count * 100
#   }
#
#   return(all_lins)
#
# }

# prot_da_lin <- reactive({switch(input$linSelec,
#                                 "PspA"= pspa_DA_Lin,
#                                 "PspB"= pspb_DA_Lin,
#                                 "PspC"= pspc_DA_Lin,
#                                 "LiaG"= liag_DA_Lin,
#                                 "LiaF"= liaf_DA_Lin,
#                                 "LiaI"= liai_DA_Lin
# )})
# prot<- reactive({
#   prot_da_lin%>% group_by(DomArch.norep)%>%
#     summarise(totalcount = sum(count)) %>%
#     filter(totalcount >1) %>% arrange(desc(totalcount))
# })
#
# prot_filt <- reactive({
#   filter(prot,totalcount >= input$cutoff)
# })
# current_row <- reactive({
#   length(prot_filt$totalcount)
# }
# )
# current_cutoff <- reactive({
#   prot[current_row,"totalcount"]
# })
#   observeEvent(input$nr,{
#     updateSliderInput(session,inputId = "cutoff",value = current_cutoff)
#   })
#
#   observeEvent(input$cutoff,  {
#     updateSliderInput(session = session, inputId = "nr", value = current_row)
#   })


#   observe({
#     current_row <- as.numeric(input$nr)
#     prot_da_lin <- switch(input$linSelec,
#            "PspA"= pspa_DA_Lin,
#            "PspB"= pspb_DA_Lin,
#            "PspC"= pspc_DA_Lin,
#            "LiaG"= liag_DA_Lin,
#            "LiaF"= liaf_DA_Lin,
#            "LiaI"= liai_DA_Lin
#            )
#
#
#     prot <- prot_da_lin %>% group_by(DomArch.norep)%>%
#       summarise(totalcount = sum(count)) %>%
#       filter(totalcount >1) %>% arrange(desc(totalcount))
#     cutoff_val <- prot[current_row,"totalcount"]
#   })
#
#   prot<- reactive({
#     prot_da_lin <- switch(input$linSelec,
#                           "PspA"= pspa_DA_Lin,
#                           "PspB"= pspb_DA_Lin,
#                           "PspC"= pspc_DA_Lin,
#                           "LiaG"= liag_DA_Lin,
#                           "LiaF"= liaf_DA_Lin,
#                           "LiaI"= liai_DA_Lin
#     )
#     prot <- prot_da_lin%>% group_by(DomArch.norep)%>%
#       summarise(totalcount = sum(count)) %>%
#       filter(totalcount >1) %>% arrange(desc(totalcount))
#   })
#   cutoff_val <- reactive({
#     prot[current_row,"totalcount"]}
#   )
#
#   observeEvent(input$cutoff, {
#     v$x <- input$cutoff
#   })
#   observeEvent(input$nr, {
#     v$x <- prot[input$nr,"totalcount"]
#   })
#
# observeEvent(input$cutoff,{
#   if(v$x != input$cutoff) updateSliderInput(session, "cutoff",value=cutoff_val)
#
#   if(v$x != input$nr) updateSliderInput(session,"nr",value =length(filter(prot,totalcount>=cutoff_val)))
#   }
#)