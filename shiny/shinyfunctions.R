library(tidyverse)
library(UpSetR)
library(gridExtra)
library(docstring)


#input a table of type xxx_GC_Lin for lineagecounts
#input a cutoff value for cutoff. A table with all total counts of values greater than cutoff will be generated

replace_doms_all <- function(prot,domains_rename, domains_remove){
  DomArch.old <- prot$DomArch

  #replace domains based on the domains_rename list
  for(j in 1:length(domains_rename$old)){
    prot <- map(prot,function(x) x %>% str_replace_all(as.vector(domains_rename$old[j]),as.vector(domains_rename$new[j])))
  }
  # #remove domains based on the domains_remove list
  for(j in 1:length(as.vector(domains_remove$domains))){
    prot <- map(prot,function(x)x %>%str_remove_all(as.vector(domains_remove$domains[j])))
  }
  #remove '+' at the start and end, as well as consecuative '+'
  prot <- map(prot, function(x) x%>%  str_replace_all("\\++\\+","\\+")
              %>% str_replace_all("^\\+","") %>%
                str_replace_all("\\+$",""))

                       # str_replace_all("\\+", " ") %>%
                       # str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
                       # str_replace_all(" ", "+")) %>% as.data.frame()

  return(as.data.frame(prot))
}

rd <- replace_doms_all(all_op_ins,domains.replace,domains.remove)

#Replace twos for all dataframe
# ab<-map(all_op_ins,function(x) x %>%  		str_replace_all("\\+", " ") %>%
#          str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
#          str_replace_all(" ", "+")) %>% as.data.frame()

#repeats2s <- function(x){
  # 	x %>%
  		# str_replace_all("\\+", " ") %>%
  		# str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
  		# str_replace_all(" ", "+")
  # }

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


#very similar to upset.plot but with the wc param. This param is used for the word count file
lineage.upset <- function(query_data,cutoff,type, wc){

  switch(type,
         da2doms=colname <- "DomArch.norep",
         gc2da= colname <- "GenContext.norep")

  words.gecutoff <- filter(wc, freq>=cutoff)

  ## Create columns for domains/DAs and fill them with 1/0
  for(i in words.gecutoff$words){
    j <-str_replace_all(string=i, pattern="\\(",replacement="\\\\(")
    j <- str_replace_all(string=j, pattern="\\)", replacement="\\\\)")
    j <- str_replace_all(string=j, pattern="\\+", replacement="\\\\+")
    j <- str_replace_all(string=j, pattern="\\_",replacement="\\\\_")
    query_data[[i]]<-if_else(grepl(j,as.matrix(query_data[,colname])),
                              true=1, false=0)
  }

  ## Creating UpSet data
  upset <- query_data %>%
    select(AccNum, Lineage, GenContext.norep, DomArch.norep,words.gecutoff$words) %>%
    mutate_all(list( ~ if(is.numeric(.))as.integer(.) else .))  %>%
    as.data.frame()


  ## Fix order of x and y variables
  upset.cutoff <-upset %>% within(colname <- factor(colname, levels=names(sort(table(colname),
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






# #########
# ##MSA
# #########
# #This will create the Phylogenetic tree
# #Pass the filepath of the fastafile as well as one of three plot types:
# #1.apeTree , 2.ggTree, 3. msaTree
# phylo.plots <- function(filepath, plot_type){
#   fasta_filepath <- "data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt"
#   my_seqs <- readAAStringSet(fasta_filepath)
#   my_seqs_msa <- msa(my_seqs)
#
#   my_seqs_msa_aln <- msaConvert(my_seqs_msa, type="seqinr::alignment")
#
#   d <- dist.alignment(my_seqs_msa_aln, "identity")
#   #  as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]
#
#   ## Phylogenetic tree
#   ## using package ape
#   ## build neighbor-joining tree
#   seqTree <- nj(d)
#   if(plot_type== "apeTree"){
#     plot(seqTree, main="Phylogenetic Tree of MSA")
#   }
#   else{
#     ## drawing trees using ggtree
#     #ggtree(seqTree)
#     groupInfo <- split(seqTree$tip.label,
#                        gsub("_\\w+", "", seqTree$tip.label))
#     seqTree <- groupOTU(seqTree, groupInfo)
#     if(plot_type == "ggTree"){
#       ggtree(seqTree, aes(color=group),
#              layout='circular') +
#         geom_tiplab(size=1, aes(angle=angle))
#     }
#     else if(plot_type == "msaTree"){
#       ## Tree + MSA
#       # beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
#       # beast_tree <- read.beast(beast_file)
#       # fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
#       # msaplot(ggtree(beast_tree), fasta)
#       offs <- 0
#       msaplot(ggtree(seqTree), fasta=fasta_filepath, offset=0.75) +
#         geom_tiplab(size=2, align=TRUE, linesize=.5)
#
#     }
#   }
# }
