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
library(docstring)

#################
## UpSet Plots
#################
upset.plot <- function(query_data="toast_rack.sub",
                       cutoff=10, type="da2doms") {
  #' UpSet Plot
  #' @author Janani Ravi
  #' @keywords UpSetR, Domains, Domain Architectures, GenomicContexts
  #' @description UpSet plot for Domain Architectures vs Domains and
  #' Genomic Contexts vs Domain Architectures.
  #' @param query_data Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is toast_rack.sub
  #' @param cutoff Numeric. Cutoff for word frequency. Default is 10.
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

  switch(type, # DA.doms.wc;
         da2doms={wc <- query_data %>% elements2words(column = "DomArch",conversion_type = type)%>%
           words2wc();
                  colname <- "DomArch"},
         gc2da={wc <- query_data %>% elements2words(column = "GenContext",conversion_type = type) %>% words2wc();
                colname <- "GenContext"})

  ## Cutoff for most/least prevalent words
  words.gecutoff <- filter(wc, freq>=cutoff) # words.gecutoff <- DA.doms.wc
  # words.ltcutoff <- filter(wc, freq<cutoff)

  ## Create columns for domains/DAs and fill them with 1/0
  for(i in words.gecutoff$words)
  {
    j <- str_replace_all(string=i, pattern="\\(", replacement="\\\\(")
    j <- str_replace_all(string=j, pattern="\\)", replacement="\\\\)")
    j <- str_replace_all(string=j, pattern="\\+", replacement="\\\\+")
    j <- str_replace_all(string=j, pattern="\\_", replacement="\\\\_")
    query_data[[i]] <- if_else(grepl(j, as.matrix(query_data[,colname])),
                              true=1, false=0)
  }
  ## Creating UpSet data
  upset <- query_data %>%
    # filter(grepl(queryname, Query)) %>%
    select(AccNum, Lineage, GenContext, DomArch,
           words.gecutoff$words) %>% ###What is this words.gecutoff$words?
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
        sets=words.gecutoff$words, sets.bar.color="turquoise3",
        main.bar.color="coral3",									#56B4E9 lightblue
        group.by="degree", order.by=c("freq"),		# "degree"
        mb.ratio=c(0.3, 0.7), # nintersects=20,
        number.angles=0, point.size=2, line.size=0.8,
        show.numbers="yes", shade.alpha = 0.25,
        mainbar.y.label="Intersection counts",
        sets.x.label="Individual set counts",
        query.legend="top")
}

###################
## Lineage Plots ##
###################

lineage.DA.plot <- function(query_data="prot",
                            #query.summ.byLin="prot.DA.summ.byLin",
                            colname="DomArch",
                            type="da2doms", cutoff = 0){ # query.elements, query.words,
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
  #' @examples lineage.DA.plot(toast_rack_data, 10, "DomArch.norep", "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query_data$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query_data$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  switch(type,
         da2doms={colname <- "DomArch";
         query.summ.byLin <- summ.DA.byLin(query_data) %>% total_counts(cutoff = cutoff, type="DA")},
         gc2da={colname <- "GenContext";
         query.summ.byLin <- summ.GC.byLin(query_data) %>% total_counts(cutoff = cutoff, type="GC")
         })


  query_data <- query_data %>% filter(grepl("a", Lineage))

  query.summ.byLin.ggplot <- drop_na(query.summ.byLin) %>%
    filter(count>1) %>%  # count or total count?
    within(Lineage <- factor(Lineage,
                                   levels=names(sort(table(Lineage),
                                                     decreasing=TRUE)))) %>%
    within(colname <- factor(colname,
                             levels=names(sort(table(colname),
                                               decreasing=F))))
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
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
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

wordcloud_element <- function(type="da2doms",
                                # query_data="prot",
                                min_freq=10){
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

  switch(type, # DA.doms.wc;
         da2doms={
           wc <- query_data %>% elements2words(column = "DomArch",conversion_type = type) %>% words2wc();
           # colname <- "DomArch.norep"
           },
         gc2da={
           wc <- query_data %>% elements2words(column = "GenContext",conversion_type = type) %>% words2wc();
           # colname <- "GenContext.norep"
         })

  ## commented out since we are directly reading in the wordcount dataframe
  # temp <- query_data$colname %>%
  #   words2wc()

  wordcloud(wc$words, wc$freq, min.freq = min_freq,
            colors=brewer.pal(8, "Spectral"),scale=c(2.5,.5))
  # wordcloud(GC.DA.wc$words,GC.DA.wc$freq,min.freq = min_freq,
  #           colors=brewer.pal(8, "Spectral"), scale=c(2.5,.4))
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
