library(tidyverse)


clean_clust_file <- function(path, writepath=NULL,query){
  colnames_prot <- c("AccNum", "GenContext","PFAM", "DomArch","arch.TMSIG","Length" ,"GenName","Lineage","Species", "Annotation","GI")

  prot <- read_tsv(path,col_names = F)

  clust <- prot %>% filter(grepl("^#",X1))

  prot <- prot %>%
    separate(X1,colnames_prot,sep=("  +"))%>% mutate(ClustName = " ",ClustID="")

  #find Which columns have # in them
  #Add Co
  ind_with_num <- which(grepl("^#",prot$AccNum))
  clsid <- separate(prot[ind_with_num,"AccNum"],col=AccNum, into=c("ClustID","ClusName"),sep = "; ")

  for(x in 1:length(clsid$ClustID)){
    lng <- str_length(x)
    clsid$ClustID[x] <- gsub(pattern="# ", replacement = paste0(strrep("0",(6-lng)),x,"."),clsid$ClustID[x])
  }


  ind_with_num <- ind_with_num %>% append(length(prot))
  #CLS id and CLS name

  for (x in length(ind_with_num):2) {
    prot[which(as.numeric(rownames(prot))
               < ind_with_num[x]),"ClustName"]= clsid[x-1,"ClusName"]
    prot[which(as.numeric(rownames(prot))
               < ind_with_num[x]),"ClustID"]= clsid[x-1,"ClustID"]
  }

  prot <- filter(prot,!grepl("^#",prot$AccNum))%>%
    mutate(Query=query)

  prot <- select(prot,AccNum,Query,ClustID,ClustName,GenContext,PFAM,DomArch,arch.TMSIG,Length,GenName,Lineage,Species,Annotation,GI)

  write_tsv(prot, writepath)

  return(prot)
  }
