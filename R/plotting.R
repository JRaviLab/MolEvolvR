#### Functions to create UpSet & Lineage plots ####
## For Domains vs Domain Architectures
## For DomArchs vs Genomic Contexts

## Modified: Apr 25, 2019
## Created: Aug 21, 2017
## Author: Janani Ravi (@jananiravi)

#################
## Pkgs needed ##
#################
library(tidyverse)
library(UpSetR)
library(gridExtra)
library(wordcloud)
library(wordcloud2)
library(sunburstR)
library(d3r)
library(viridis)

shorten_lineage <- function(data, colname = "Lineage")
{
  abbrv <- function(x){
    pos_gt = str_locate(x,">")
    pos_gt = pos_gt[1]
    if(is.na(pos_gt)){
      return(toupper(substr(x,1,1)))
    }
    return(paste0(toupper(substr(x,1,1)),substr(x,pos_gt, nchar(as.character(x)))))
  }
  # Shorten lineages to include only the first letter of kingdom
  data$Lineage <- unlist((pmap(list(data$Lineage), function(x) abbrv(x) )))
  return(data)
}


#################
## UpSet Plots
#################
upset.plot <- function(query_data="toast_rack.sub",
                       colname = "DomArch", cutoff = 90,
                       RowsCutoff = FALSE) {
  #' UpSet Plot
  #' @author Janani Ravi
  #' @keywords UpSetR, Domains, Domain Architectures, GenomicContexts
  #' @description UpSet plot for Domain Architectures vs Domains and
  #' Genomic Contexts vs Domain Architectures.
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is toast_rack.sub
  #' @param cutoff Numeric. Cutoff for word frequency. Default is 90.
  #' @param type Character. Either "da2doms" for Domains vs Domain Architectures
  #' or "gc2da" for Domain Architectures (of neighbors) vs Genomic Contexts.
  #' Default is "da2doms"
  #' @examples upset.plot(pspa.sub, 10, "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  # Get Total Counts
  # colname = string(colname)
  tc <- query_data %>% total_counts(column =  colname, cutoff = cutoff, RowsCutoff = RowsCutoff, digits = 5)


  ##### Remove Tails ####
  # tails comprise of less than 1% of data each
  # ie) individual percent is less than 1
  # Note: this is based on frequency of the actual DA's and GC's, not on the frequency of the constituent words


  #####


  column <- sym(colname)

  if(grepl("GenContext", colname, ignore.case = T))
  {
    type = "gc2da"
  }
  else if(grepl("DomArch|ClustName", colname, ignore.case = T))
  {
    type = "da2doms"
  }

  # Get words from filter
  words.tc <- tc %>% select({{column}}) %>% distinct() %>%
    elements2words(column = colname,conversion_type = type)
  # names(words.tc)[1] <- "words"
  words.tc <- words.tc %>% str_split(pattern = " ")
  words.tc = as.data.frame(words.tc, col.names = "Words", stringsAsFactors = F) %>% filter(!grepl("^X$|^X\\(s\\)$",Words)) %>% pull(Words)# remove "X" and "X(s)"
  words.tc <- words.tc[2:(length(words.tc)-1)]

  # Only use the DAs/GCs that are within the total count cutoff
  query_data <- query_data %>% filter({{column}} %in% (select(tc, {{column}}) %>% distinct() %>% pull({{column}}) ) )

  ## Create columns for domains/DAs and fill them with 1/0
  for(i in words.tc) #$words)
  {
    j <- str_replace_all(string=i, pattern="\\(", replacement="\\\\(")
    j <- str_replace_all(string=j, pattern="\\)", replacement="\\\\)")
    j <- str_replace_all(string=j, pattern="\\+", replacement="\\\\+")
    j <- str_replace_all(string=j, pattern="\\_", replacement="\\\\_")
    j <- str_replace_all(string=j, pattern="\\?", replacement="\\\\?")

    # Don't pick up anything with a '(' immediatly following j
    # Makes sure PspA doesn't pick up PspA(s)
    if(grepl("DomArch|ClustName", colname, ignore.case = T))
    {
      j <- paste0(j,"(?!\\()")
    }

    else if(grepl("GenContext", colname, ignore.case = T))
    {

      # Negative lookahead and lookbehind for +
      j <- paste0("(?<!\\+)",j, "(?!\\+)")
    }

    # Grep doesn't work with string+negative lookahead?
    # Use str_detect instead
    # query_data[[i]] <- if_else(grepl(j, as.matrix(query_data[,colname])),
    #                            true=1, false=0)

    query_data[[i]] <- if_else(str_detect(string = as.matrix(query_data[,colname]), pattern = j),
                         true=1, false=0)
  }

  ## Creating UpSet data
  upset <- query_data %>%
    select(AccNum, Lineage, {{column}}, # GenContext, DomArch,
           words.tc) %>%
    mutate_all(list(~if(is.numeric(.)) as.integer(.) else .)) %>%
    as.data.frame()



  ## Fix order of x and y variables
  upset.cutoff <- upset %>%
    # filter(!grepl(paste(words.ltcutoff$words, collapse="|"), colname)) %>%
    within(colname <- factor(colname, levels=names(sort(table(colname),
                                                        decreasing=TRUE))))
  ## UpSetR plot
  par(oma=c(5,5,5,5), mar=c(3,3,3,3))
  upset(upset.cutoff[c(3,5:ncol(upset.cutoff))],	# text.scale=1.5,
        #sets=words.tc,
        nsets = length(words.tc),
        nintersects = NA,
        sets.bar.color="turquoise3",
        main.bar.color="coral3",									#56B4E9 lightblue
        group.by="degree", order.by=c("freq"),		# "degree"
        mb.ratio=c(0.3, 0.7), # nintersects=20,
        number.angles=0, point.size=2, line.size=0.8,
        show.numbers="yes", shade.alpha = 0.25,
        mainbar.y.label="Intersection counts",
        sets.x.label="Individual set counts",
        query.legend="top",
        text.scale = 1.5
        )
}

###################
## Lineage Plots ##
###################

lineage.DA.plot <- function(query_data="prot",
                            colname="DomArch",
                            cutoff = 90,
                            RowsCutoff = FALSE,
                            color = "default"){ # query.elements, query.words,
  #' Lineage Plot: Heatmap of Domains/DAs/GCs vs Lineages
  #' @author Janani Ravi
  #' @keywords Lineages, Domains, Domain Architectures, GenomicContexts
  #' @description Lineage plot for Domains, Domain Architectures and
  #' Genomic Contexts. Heatmap.
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format).
  #' Default is prot (variable w/ protein data).
  #' @param query.summ.byLin Output of summ.DA.byLin(XXX.sub) or summ.GC.byLin(XXX.sub)
  #' @param colname Column name from query_data: "DomArch.norep", "GenContext.norep",
  #' "DomArch.PFAM.norep" or "DomArch.LADB.norep". Default is "DomArch.norep".
  #' @param type Character. Default is "da2doms" for Domain Architectures.
  #' Other alternative: "gc2da" for Genomic Contexts.
  #' @param color Color for the heatmap. One of six options: "default", "magma", "inferno",
  #' "plasma", "viridis", or "cividis"
  #' @examples lineage.DA.plot(toast_rack_data, 10, "DomArch.norep", "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.


  query_data <- shorten_lineage(query_data, "Lineage")

  query.summ.byLin <- query_data %>% total_counts(cutoff = cutoff, column = colname, RowsCutoff = RowsCutoff)

  query.summ.byLin$Lineage <- map(query.summ.byLin$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
    unlist()

  query_data <- query_data %>% filter(grepl("a", Lineage))

  column <- sym(colname)

  query.summ.byLin.ggplot <-  drop_na(query.summ.byLin) %>%
                                filter(count>1) %>%# count or total count?
                                arrange(totalcount)

  query.summ.byLin.ggplot$Lineage <- factor(query.summ.byLin.ggplot$Lineage,
                                            levels = sort(names(sort(table(query.summ.byLin.ggplot$Lineage),decreasing=TRUE))))

  query.summ.byLin.ggplot[,colname] <- factor(pull(query.summ.byLin.ggplot, {{column}}) ,
                                            levels = (pull(query.summ.byLin.ggplot, {{column}}) %>% unique()))


  if(color == "default"){
  ## Tile plot
  ggplot(data=query.summ.byLin.ggplot,
         aes_string(x="Lineage", y=colname)) +
    geom_tile(data=subset(query.summ.byLin.ggplot,
                          !is.na(count)),
              aes(fill=count),
              colour="darkred", size=0.3) + #, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    scale_x_discrete(position="top") +
    theme_minimal() + # coord_flip() +
    theme(axis.text.x=element_text(angle=65,hjust=0,vjust=0.5))
  }
  else
  {
    ggplot(data=query.summ.byLin.ggplot,
           aes_string(x="Lineage", y=colname)) +
      geom_tile(data=subset(query.summ.byLin.ggplot,
                            !is.na(count)),
                aes(fill=count))+
                #colour="darkred", size=0.3) + #, width=0.7, height=0.7),
      scale_fill_viridis(discrete = F, option = color)+
      scale_x_discrete(position="top") +
      theme_minimal() + # coord_flip() +
      theme(axis.text.x=element_text(angle=65,hjust=0,vjust=0.5))
  }
}



lineage.Query.plot <- function(query_data=all,
                               queries,
                               colname = "ClustName",
                               cutoff){
  #' Lineage Plot: Heatmap of Queries vs Lineages
  #' @authors Janani Ravi, Samuel Chen
  #' @keywords Lineages, Domains, Domain Architectures, GenomicContexts
  #' @description Lineage plot for queries. Heatmap.
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format).
  #' Default is prot (variable w/ protein data).
  #' @param queries Character Vector containing the queries that will be used for the categories
  #' @examples lineage.Query.plot(prot, c("PspA","PspB","PspC","PspM","PspN"), 95)
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  query_by_lineage <- function(data, query, column, by){
    column <- sym(column); by <- sym(by)

    # filter the protein by the query
    data <- data %>% filter(grepl(pattern=query, x={{column}},
                                  ignore.case=T)) %>% select({{by}},count)

    if(nrow(data) == 0)
    {
      data$Query <- character(0)
      #data$count <- integer(0)
      return(data)
    }

    data$Query <- query

    data <- data %>% filter(!grepl("^-$", {{by}})) %>%
      group_by(Query, {{by}}) %>%
      summarise(count=sum(count)) %>%
      arrange(desc(count))
    return(data)

  }

  col <- sym(colname)

  query_data <- query_data %>% total_counts(column = colname, cutoff = cutoff)

  # query_data contains all rows that possess a lineage
  query_data <- query_data %>% filter(grepl("a", Lineage))

  query_data <- shorten_lineage(query_data, "Lineage")

  query_lin_counts = data.frame("Query" = character(0), "Lineage" = character(0), "count"= integer())
  for(q in queries){
    query_lin <- query_by_lineage(data = query_data, query = q, column = {{col}}, by = "Lineage")
    query_lin_counts <- dplyr::union(query_lin_counts, query_lin)
  }

  query_lin_counts$Lineage <- map(query_lin_counts$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
    unlist()

  query.summ.byLin.ggplot <- drop_na(query_lin_counts) %>%
    filter(count>1) %>%  # count or total count?
    within(Lineage <- factor(Lineage,
                             levels= sort(names(sort(table(Lineage),
                                               decreasing=TRUE))))) %>%
    within(Query <- factor(Query,
                           levels=names(sort(table(Query),
                                             decreasing=F))))
  ## Tile plot
  ggplot(data=query.summ.byLin.ggplot,
         aes_string(x="Lineage", y="Query")) +
    geom_tile(data=subset(query.summ.byLin.ggplot,
                          !is.na(count)),
              aes(fill=count),
              colour="darkred", size=0.3) + #, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    scale_x_discrete(position="top") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=65,hjust=0,vjust=0.5))
  # scale_y_discrete(position = "bottom") +
  # scale_x_discrete(position = "bottom") +
  # coord_flip()
}


lineage.neighbors.plot <- function(query_data="prot", query="pspa",
                                   colname="GenContext.norep"){
  #' Lineage Plot for top neighbors
  #' @author Janani Ravi
  #' @keywords Lineages, Domains, Domain Architectures, GenomicContexts
  #' @description Lineage plot for top neighbors obtained from DAs of
  #' Genomic Contexts.
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is pspa_data.
  #' @param query Name of query protein/domain. Default is "pspa".
  #' @param colname Column name from query_data. Default is "GenContext.norep".
  #' @examples lineage.neighbors.plot(pspa_data, pspa, "GenContext.norep", "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  query_data <- query_data %>% filter(grepl("a", Lineage))
  query.GCDA <- read_delim(paste0("Top-",query,"-neighbors.txt"),
                           delim="\t", escape_double=FALSE, col_names=FALSE,
                           na="NA", comment="#", trim_ws=TRUE)
  colnames(query.GCDA) <- "domarchs"

  neighbors.da <- query.GCDA$domarchs
  for(i in neighbors.da)
  {
    j <- str_replace_all(string=i, pattern="\\(", replacement="\\\\(")
    j <- str_replace_all(string=j, pattern="\\)", replacement="\\\\)")
    j <- str_replace_all(string=j, pattern="\\+", replacement="\\\\+")
    j <- str_replace_all(string=j, pattern="\\_", replacement="\\\\_")
    # j <- str_replace_all(string=j, pattern="\\-", replacement="\\\\-")
    query_data[[i]] <- if_else(grepl(j, as.matrix(query_data[,colname])),
                               true=1, false=0)
  }

  query.ggplot <- query_data %>%
    gather(key=TopNeighbors.DA, value=count, 19:ncol(query_data)) %>%
    select("Lineage", "TopNeighbors.DA", "count") %>% # "DomArch.norep","GenContext.norep",
    group_by(TopNeighbors.DA, Lineage) %>%
    summarise(lincount=sum(count), bin=as.numeric(as.logical(lincount))) %>%
    arrange(desc(lincount)) %>%
    within(TopNeighbors.DA <- factor(TopNeighbors.DA,
                                     levels=rev(names(sort(table(TopNeighbors.DA),
                                                           decreasing=TRUE))))) %>%
    within(Lineage <- factor(Lineage,
                             levels=names(sort(table(Lineage),
                                               decreasing=TRUE))))

  ggplot(query.ggplot, aes(x=Lineage, y=TopNeighbors.DA)) +
    geom_tile(data=subset(query.ggplot,
                          !is.na(lincount)),		# bin
              aes(fill=lincount),								# bin
              colour="coral3", size=0.3) +			#, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    scale_x_discrete(position="top") +				# guides(fill=FALSE) +
    theme_minimal() +														# coord_flip() +
    theme(axis.text.x=element_text(angle=90,
                                   hjust=0,vjust=0.5))
}

lineage.domain_repeats.plot <- function(query_data, colname) {
  # query_data <- pspa_data
  # colname <- "SIG.TM.LADB"

  ## Create columns for domains/DAs and fill them with 1/0
  for(i in query.DAdoms$domains)
  {
    j <- str_replace_all(string=i, pattern="\\(", replacement="\\\\(")
    j <- str_replace_all(string=j, pattern="\\)", replacement="\\\\)")
    j <- str_replace_all(string=j, pattern="\\+", replacement="\\\\+")
    j <- str_replace_all(string=j, pattern="\\_", replacement="\\\\_")
    # query_data[[i]] <- if_else(grepl(j, as.matrix(query_data[,colname])),
    # 													true=1, false=0)		## BINARY
    query_data[[i]] <- str_count(string=as.matrix(query_data[,colname]),
                                 pattern=j)					## ACTUAL COUNTS
  }

  ## Subsetting relevant columns
  ggplot.data <- query_data %>%
    # filter(grepl(queryname, Query)) %>%
    select(DomArch.norep, Lineage, GenContext.norep,
           SIG.TM.LADB, GenContext, AccNum,
           query.DAdoms$domains) %>% # words.gecutoff$words
    mutate_all(list(~ if(is.numeric(.)) as.integer(.) else .)) %>%
    as.data.frame()

  # ## written on Sep 4, 2017
  # write_delim(x=ggplot.data, "toast-rack.domain_repeat_counts.v1-2.txt",
  # 						delim="\t", col_names=TRUE)

  ## Gathering element/word columns
  ggplot.data.gather <- ggplot.data %>%
    gather(key=domains, value=count, 7:ncol(ggplot.data)) # %>%
  # select(DomArch.norep, Lineage, domains, count)

  # ## written on Sep 4
  # write_delim(ggplot.data.gather,
  # 						"toast-rack.domain_repeat_counts-gathered.v1-2.txt",
  # 						delim="\t", col_names=TRUE)

  ## Stacked column plot
  ggplot(data=ggplot.data.gather, aes(x=Lineage, y=domains)) + # aes_string # plot <- (
    # geom_col(position="fill") +
    geom_tile(data=subset(ggplot.data.gather, !is.na(count)),
              aes(fill=as.numeric(as.logical(count))),
              colour="coral3", size=0.3) + #, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    scale_x_discrete(position="top") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
}


################
## Wordclouds ##
################
#### NEEDS SOME WORK


# Modified script of wordcloud2
wordcloud3 <- function (data, size = 1, minSize = 0, gridSize = 0, fontFamily = "Segoe UI",
                        fontWeight = "bold", color = "random-dark", backgroundColor = "white",
                        minRotation = -pi/4, maxRotation = pi/4, shuffle = TRUE,
                        rotateRatio = 0.4, shape = "circle", ellipticity = 0.65,
                        widgetsize = NULL, figPath = NULL, hoverFunction = NULL) {
  if ("table" %in% class(data)) {
    dataOut = data.frame(name = names(data), freq = as.vector(data))
  }
  else {
    data = as.data.frame(data)
    dataOut = data[, 1:3]
    names(dataOut) = c("name", "freq", "label")
  }
  if (!is.null(figPath)) {
    if (!file.exists(figPath)) {
      stop("cannot find fig in the figPath")
    }
    spPath = strsplit(figPath, "\\.")[[1]]
    len = length(spPath)
    figClass = spPath[len]
    if (!figClass %in% c("jpeg", "jpg", "png", "bmp", "gif")) {
      stop("file should be a jpeg, jpg, png, bmp or gif file!")
    }
    base64 = base64enc::base64encode(figPath)
    base64 = paste0("data:image/", figClass, ";base64,",
                    base64)
  }
  else {
    base64 = NULL
  }
  weightFactor = size * 180/max(dataOut$freq)
  settings <- list(word = dataOut$name, freq = dataOut$freq, label = dataOut$label,
                   fontFamily = fontFamily, fontWeight = fontWeight, color = color,
                   minSize = minSize, weightFactor = weightFactor, backgroundColor = backgroundColor,
                   gridSize = gridSize, minRotation = minRotation, maxRotation = maxRotation,
                   shuffle = shuffle, rotateRatio = rotateRatio, shape = shape,
                   ellipticity = ellipticity, figBase64 = base64, hover = htmlwidgets::JS(hoverFunction))
  chart = htmlwidgets::createWidget("wordcloud2", settings,
                                    width = widgetsize[1], height = widgetsize[2], sizingPolicy = htmlwidgets::sizingPolicy(viewer.padding = 0,
                                                                                                                            browser.padding = 0, browser.fill = TRUE))

  chart
  #htmlwidgets::onRender(chart, "function(el,x){\n                        console.log(123);\n                        if(!iii){\n                          window.location.reload();\n                          iii = False;\n\n                        }\n  }")
}


wordcloud_element <- function(query_data="prot",
                              colname = "DomArch",
                              cutoff = 70,
                              UsingRowsCutoff = FALSE
){
  #' Wordclouds for the predominant domains, domain architectures.
  #' @author Janani Ravi
  #' @keywords Domains, Domain Architectures, GenomicContexts
  #' @description Wordclouds for the predominant domains (from DAs) and DAs (from GC)
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is "prot".
  #' @param type Character. Default is "da2doms" for Domain Architectures.
  #' Other alternative: "gc2da" for Genomic Contexts.
  #' @examples wordcloud_element(prot, "da2doms", 10)
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  # Get Total Counts
  # colname = string(colname)

  tc <- query_data %>% total_counts(column =  colname, cutoff = cutoff, RowsCutoff = UsingRowsCutoff, digits = 5)

  column <- sym(colname)
  # Get words from filter
  # words.tc <- tc %>% select({{column}}, totalcount) %>% distinct()

  query_data = query_data %>% filter({{column}} %in% pull(tc,{{colname}}))

  if(grepl("DomArch|Clustname", colname))
  {
    type = "da2doms"
  }
  else if(grepl("GenContext", colname))
  {
    type = "gc2da"
  }

  words.tc <- query_data %>% elements2words(column = colname,
                                            conversion_type = type) %>%
    words2wc()

  # names(words.tc) <- c("words", "freq")

  # need a label column for actual frequencies, and frequencies will be the
  # normalized sizes
  # words.tc$label <- words.tc$freq
  #
  # words.tc <- words.tc %>% mutate(freq = log10(freq))
  #
  # words.tc <- words.tc %>% select(words, freq, label)
  #
  # wordcloud3(words.tc, minSize = 0)
  wordcloud(words.tc$words, words.tc$freq, min.freq = 1,
            colors=brewer.pal(8, "Spectral"),scale=c(4.5,1))
}


wordcloud2_element <- function(query_data="prot",
                              colname = "DomArch",
                              cutoff = 70,
                              UsingRowsCutoff = FALSE
){
  #' Wordclouds for the predominant domains, domain architectures.
  #' @author Janani Ravi
  #' @keywords Domains, Domain Architectures, GenomicContexts
  #' @description Wordclouds for the predominant domains (from DAs) and DAs (from GC)
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is "prot".
  #' @param type Character. Default is "da2doms" for Domain Architectures.
  #' Other alternative: "gc2da" for Genomic Contexts.
  #' @examples wordcloud_element(prot, "da2doms", 10)
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.
  tc <- query_data %>% total_counts(column =  colname, cutoff = cutoff, RowsCutoff = UsingRowsCutoff, digits = 5)

  column <- sym(colname)
  query_data = query_data %>% filter({{column}} %in% pull(tc,{{colname}}))

  if(grepl("DomArch|Clustname", colname))
  {
    type = "da2doms"
  }
  else if(grepl("GenContext", colname))
  {
    type = "gc2da"
  }

  words.tc <- query_data %>% elements2words(column = colname,
                                            conversion_type = type) %>%
    words2wc()

  names(words.tc) <- c("words", "freq")

  # need a label column for actual frequencies, and frequencies will be the
  # normalized sizes
  words.tc$label <- words.tc$freq

  words.tc <- words.tc %>% mutate(freq = log10(freq))

  words.tc <- words.tc %>% select(words, freq, label)

  wordcloud3(words.tc, minSize = 0)
}





#### Sunburst #####
lineage_sunburst <- function(prot, lineage_column = "Lineage",
                             type = "sunburst",
                             levels = 2, colors = NULL, legendOrder = NULL)
{
  #'
  #'
  #'@param prot Data frame containing a lineage column that the sunburst plot will be generated for
  #'@param lineage_column String. Name of the lineage column within the data frame. Defaults to "Lineage"
  #'@param type String, either "sunburst" or "sund2b". If type is "sunburst", a sunburst plot of the lineage
  #'@param levels Integer. Number of levels the sunburst will have.
  #'@param legendOrder String vector. The order of the legend. If legendOrder is NULL,
  #' then the legend will be in the descending order of the top level hierarchy.
  #'will be rendered. If the type is sund2b, a sund2b plot will be rendered.


  lin_col <- sym(lineage_column)

  levels_vec = c()
  for(i in 1:levels)
  {
    levels_vec = append(levels_vec, paste0("level",i))
  }

  # Take lineage column and break into the first to levels
  prot <- prot %>% select({{lin_col}})
  protLevels <- prot %>% separate({{lin_col}}, into = levels_vec, sep = ">")
  # Count the occurrance of each group of levels
  protLevels = protLevels %>% group_by_at(levels_vec) %>% summarise(size = n())

  tree <- d3_nest(protLevels, value_cols = "size")

  # Plot sunburst
  if(type == "sunburst")
  {
    sunburst(tree, legend = list(w = 150, h = 15, r = 3, s = 3), colors = colors, legendOrder = legendOrder)
  }
  else if(type == "sund2b")
  {
    sund2b(tree)
  }
}


## COMMENTED LINEAGE.DA.PLOT
# lineage.plot <- function(query_data, cutoff, type) {
# 	switch(type,
# 				 da2doms={wc <- DA.doms.wc; words <- toast_rack.DAdoms; colname <- "DomArch.norep"; toast_rack.summ.byLin <- toast_rack.DA.summ.byLin}, # elements <- toast_rack.DA;
# 				 gc2da={wc <- GC.DA.wc; words <- toast_rack.GCDA; colname <- "GenContext.norep"; toast_rack.summ.byLin <- toast_rack.GC.summ.byLin} # elements <- toast_rack.GC;
# 	)
#
# 	# # ## Cutoff for most/least prevalent words (domains for DAs; DAs for GCs)
# 	# # words.gecutoff <- filter(wc, freq>=cutoff) # NOT NEEDED
# 	# # words.ltcutoff <- filter(wc, freq<cutoff)
# 	#
# 	# ## Create columns for domains/DAs and fill them with 1/0
# 	# # for(i in words.gecutoff$words)
# 	# for(i in elements$domarchs)
# 	# {
# 	# 	j <- str_replace(string=i, pattern="\\(", replacement="\\\\(")
# 	# 	j <- str_replace(string=j, pattern="\\)", replacement="\\\\)")
# 	# 	j <- str_replace(string=j, pattern="\\+", replacement="\\\\+")
# 	# 	j <- str_replace(string=j, pattern="\\_", replacement="\\\\_")
# 	# 	# query_data[[i]] <- if_else(grepl(j, as.matrix(query_data[,colname])),
# 	# 	# 													true=1, false=0)		## BINARY
# 	# 	query_data[[i]] <- str_count(string=as.matrix(query_data[,colname]),
# 	# 															pattern=j)					## ACTUAL COUNTS
# 	# }
# 	#
# 	# ## Subsetting relevant columns
# 	# ggplot.data <- query_data %>%
# 	# 	# filter(grepl(queryname, Query)) %>%
# 	# 	select(AccNum, Lineage.final, GenContext.norep, DomArch.norep,
# 	# 				 words) %>% # words.gecutoff$words
# 	# 	mutate_all(funs(if(is.numeric(.)) as.integer(.) else .)) %>%
# 	# 	as.data.frame()
# 	# ## Gathering element/word columns
# 	# ggplot.data.gather <- ggplot.data %>%
# 	# 	gather(key=domains, value=counts, 5:ncol(ggplot.data))
# 	# ## Summarizing data
# 	# ggplot.data.summ <- ggplot.data %>%
# 	# 	group_by(DomArch.norep, Lineage.final) %>%
# 	# 	summarise(counts=n(), counts2=n_distinct(Lineage.final)) %>%
# 	# 	filter(counts>10)
# toast_rack.summ.byLin.ggplot <- drop_na(toast_rack.summ.byLin) %>%
# 	filter(count>1) %>%
# 	within(Lineage.final <- factor(Lineage.final,
# 													 levels=names(sort(table(Lineage.final),
# 													 									decreasing=TRUE)))) %>%
# 	within(DomArch.norep <- factor(DomArch.norep,
# 																 levels=names(sort(table(DomArch.norep),
# 																 									decreasing=F))))
# 	## Tile plot
# 	## data=ggplot.data.summ
# 	ggplot(data=toast_rack.DA.summ.byLin.ggplot,
# 				 aes(x=Lineage.final, y=DomArch.norep)) + # aes_string # plot <- (
# 		geom_tile(data=subset(toast_rack.DA.summ.byLin.ggplot,
# 													!is.na(count)),
# 							aes(fill=count),
# 							colour="darkred", size=0.3) + #, width=0.7, height=0.7),
# 		scale_fill_gradient(low="white", high="darkred") +
# 		scale_x_discrete(position="top") +
# 		theme_minimal() + # coord_flip() +
# 		theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
#
# 	# ## Stacked column plot
# 	# ggplot(data=ggplot.data.gather, aes(x=domains, y=counts, fill=Lineage.final)) + # aes_string # plot <- (
# 	# 	geom_col(position="fill") +
# 	# 	# geom_tile(data=subset(query.list.lin.summ.ge5, !is.na(count)),
# 	# 						# aes(fill=counts),
# 	# 						# colour="darkred", size=0.3) + #, width=0.7, height=0.7),
# 	# 	# scale_fill_gradient(low="white", high="darkred") +
# 	# 	# scale_x_discrete(position="top") +
# 	# 	theme_minimal() +
# 	# 	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
# 	# 				axis.text.y=element_text(angle=90,hjust=1,vjust=0.5))
#
# }
