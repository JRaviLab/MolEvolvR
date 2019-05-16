library(tidyverse)
library(UpSetR)
library(gridExtra)
library(docstring)
library(Biostrings)
#for msa
library(msa)
# library(seqinr)
# library(phylogram)
# library(phytools) # depends on phangorn which depends on quadprog. so!!!
library(ape)
library(ggtree)
# library(tidytree)

#input a table of type xxx_GC_Lin for lineagecounts
#input a cutoff value for cutoff. A table with all total counts of values greater than cutoff will be generated
cummulative.count <- function(lineagecounts ,cutoff = 20, DA_GC = "GC"){
  if(DA_GC == "GC"){
    gc_count <- lineagecounts %>% group_by(GenContext.norep) %>% summarise(totalcount = sum(count))  %>% filter(totalcount >= cutoff)
    cummulative_count <- left_join(lineagecounts,gc_count, by = "GenContext.norep")
  }
  else{
    da_count <- lineagecounts %>% group_by(DomArch.norep) %>% summarise(totalcount = sum(count))  %>% filter(totalcount >= cutoff)
    cummulative_count <- left_join(lineagecounts,da_count, by = "DomArch.norep")
  }
  return(cummulative_count)
}

#Generates Heatmap
#Samuel Chen added cutoff parameter so the slider input would work
lineage.DA.plot <- function(query.sub="toast_rack.sub",
                            query.summ.byLin="toast_rack.DA.summ.byLin",
                            colname="DomArch.norep", type="da2doms",cutoff){ # query.elements, query.words,
  #' Lineage Plot: Heatmap of Domains/DAs/GCs vs Lineages
  #' @author Janani Ravi
  #' @keywords Lineages, Domains, Domain Architectures, GenomicContexts
  #' @description Lineage plot for Domains, Domain Architectures and
  #' Genomic Contexts. Heatmap.
  #' @param query.sub Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is toast_rack.sub.
  #' @param query.summ.byLin Output of summ.DA.byLin(XXX.sub) or summ.GC.byLin(XXX.sub)
  #' @param colname Column name from query.sub: "DomArch.norep", "GenContext.norep",
  #' "DomArch.PFAM.norep" or "DomArch.LADB.norep". Default is "DomArch.norep".
  #' @param type Character. Default is "da2doms" for Domain Architectures.
  #' Other alternative: "gc2da" for Genomic Contexts.
  #' @examples lineage.DA.plot(toast_rack.sub, 10, "DomArch.norep", "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query.sub$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query.sub$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  #Can accept both DA and GC
  switch(type,
         da2doms={colname <- "DomArch.norep"},
         gc2da={colname <- "GenContext.norep"})

  #Filter by column 'lineage.final', looking for a character "a".
  #Is this to get rid of NA rows?
  query.sub <- query.sub %>% filter(grepl("a", Lineage.final))# %>%  filter(totalcount >= cutoff)


  query.summ.byLin.ggplot <- drop_na(query.summ.byLin) %>%
    filter(totalcount >= cutoff) %>%
    filter(count>1) %>%
    within(Lineage.final <- factor(Lineage.final,
                                   levels=names(sort(table(Lineage.final),
                                                     decreasing=TRUE)))) %>%
    within(colname <- factor(colname,
                             levels=names(sort(table(colname),
                                               decreasing=F))))
  ## Tile plot
  ggplot(data=query.summ.byLin.ggplot,
         aes_string(x="Lineage.final", y=colname)) +
    geom_tile(data=subset(query.summ.byLin.ggplot,
                          !is.na(count)),
              aes(fill=count),
              colour="darkred", size=0.3) + #, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    scale_x_discrete(position="top") +
    theme_minimal() + # coord_flip() +
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
}

#################
## UpSet Plots
#################
#SZC Modify to include DA.doms.wc parameter
upset.plot <- function(query.sub="toast_rack.sub", doms.wc,
                       cutoff=10, type="da2doms") {
  #' UpSet Plot
  #' @author Janani Ravi
  #' @keywords UpSetR, Domains, Domain Architectures, GenomicContexts
  #' @description UpSet plot for Domain Architectures vs Domains and
  #' Genomic Contexts vs Domain Architectures.
  #' @param query.sub Data frame of protein homologs with the usual 11 columns +
  #' additional word columns (0/1 format). Default is toast_rack.sub
  #' @param cutoff Numeric. Cutoff for word frequency. Default is 10.
  #' @param type Character. Either "da2doms" for Domains vs Domain Architectures
  #' or "gc2da" for Domain Architectures (of neighbors) vs Genomic Contexts.
  #' Default is "da2doms"
  #' @examples upset.plot(pspa.sub, 10, "da2doms")
  #' @details For "da2doms" you would need the file DA.doms.wc as well as the
  #' column query.sub$DomArch.norep
  #'
  #' For "gc2da", you would need the file GC.DA.wc as well as the column
  #' query.sub$GenContext.norep
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.

  switch(type, # DA.doms.wc;
         da2doms={wc <- doms.wc; colname <- "DomArch.norep"},
         gc2da={wc <- doms.wc; colname <- "GenContext.norep"})

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
    query.sub[[i]] <- if_else(grepl(j, as.matrix(query.sub[,colname])),
                              true=1, false=0)
  }
  ## Creating UpSet data
  upset <- query.sub %>%
    # filter(grepl(queryname, Query)) %>%
    select(AccNum, Lineage, GenContext.norep, DomArch.norep,
           words.gecutoff$words) %>%
    mutate_all(funs(if(is.numeric(.)) as.integer(.) else .)) %>%
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
        # mb.ratio=c(0.6, 0.4), nintersects=20,
        number.angles=0, point.size=2, line.size=1,
        mainbar.y.label="Intersection counts",
        sets.x.label="Individual set counts",
        query.legend="top")
}


#########
##MSA
#########
#This will create the Phylogenetic tree
#Pass the filepath of the fastafile as well as one of three plot types:
#1.apeTree , 2.ggTree, 3. msaTree
phylo.plots <- function(filepath, plot_type){
  fasta_filepath <- "data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt"
  my_seqs <- readAAStringSet(fasta_filepath)
  my_seqs_msa <- msa(my_seqs)

  my_seqs_msa_aln <- msaConvert(my_seqs_msa, type="seqinr::alignment")

  d <- dist.alignment(my_seqs_msa_aln, "identity")
  #  as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

  ## Phylogenetic tree
  ## using package ape
  ## build neighbor-joining tree
  seqTree <- nj(d)
  if(plot_type== "apeTree"){
    plot(seqTree, main="Phylogenetic Tree of MSA")
  }
  else{
    ## drawing trees using ggtree
    #ggtree(seqTree)
    groupInfo <- split(seqTree$tip.label,
                       gsub("_\\w+", "", seqTree$tip.label))
    seqTree <- groupOTU(seqTree, groupInfo)
    if(plot_type == "ggTree"){
      ggtree(seqTree, aes(color=group),
             layout='circular') +
        geom_tiplab(size=1, aes(angle=angle))
    }
    else if(plot_type == "msaTree"){
      ## Tree + MSA
      # beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
      # beast_tree <- read.beast(beast_file)
      # fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
      # msaplot(ggtree(beast_tree), fasta)
      offs <- 0
      msaplot(ggtree(seqTree), fasta=fasta_filepath, offset=0.75) +
        geom_tiplab(size=2, align=TRUE, linesize=.5)

    }
  }
}
