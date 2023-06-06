suppressPackageStartupMessages(library(tidyverse))

#'Reverse Operons
#'
#'Reverses the direction of the operons
#'
#'This function reverses the direction of operons in a column of a data frame called 'GenContext'. Returns
#'the original dataframe with the original operons in a column called 'GenContext.old' and the reversed
#'operons in a column called 'GenContext'
#'
#'@param prot A data frame containing a GenContext column
#'@examples reverse_operons(pspa.sub)
reverse_operons <- function(prot){
  #operon straightner
  opvec <- prot$GenContext.norep   #operon column from a df containing cleanblastclust columns
  opvec <- gsub(pattern = ">",replacement = ">|",x = opvec)%>%
    gsub(pattern = "<",replacement = "|<") %>%
    gsub(pattern = "\\|\\|",replacement = "\\|=\\|")

  #op.list takes opvec and for each element in opvec, a list of characters are created by spliting the strings
  op.list <- strsplit(x = opvec, split =  "\\|")%>%
    replace_na("-")
  #Removes empty strings from the list
  op.list <- map(1:length(op.list),
                 function(x) {
                   if(any(op.list[[x]]=="")) op.list[[x]][which(op.list[[x]] !="")]
                   else op.list[[x]] })
  #te is the
  te <- map(1:length(op.list),
            function(x) op.list[[x]][grep("\\*", op.list[[x]])])
  torev <- grep("^<",te)

  te <- op.list[torev]%>%
    map( function(x) gsub(pattern = "<-|->", replacement = "",x = x)) %>%
    map(rev)
  witheq <- grep(pattern = "=",x = te)
  withouteq <- which(!((1:length(te)) %in% witheq))
  ge <- te[witheq]

  reveql=function(x){
    w <- x
    y <- rep(NA, length(w))
    d <- 1
    for(j in 1:length(w)){
      if(w[j]=="=") d <- d*(-1)
      if(d==1 && w[j] != "=") {y[j]=paste(w[j],"->", sep = "")
      } else if (d==-1 && w[j] != "="){
        y[j] <- paste("<-",w[j], sep = "")
      } else{
        y[j] <- "="
      }
    }
    return(y)
  }

  ge <- map(1:length(ge), function(x) reveql(ge[[x]]))
  ye <- te[withouteq]
  ye <- map(1:length(ye), function(x) unname(sapply(ye[[x]], function(y) paste(y,"->", sep = ""))))

  te[witheq] <- ge
  te[withouteq]<-ye
  op.list[torev] <- te

  revopvec <- unlist(map(op.list, function(x) paste(x, collapse = "")))%>%
    gsub(pattern = "=",replacement = "\\|\\|") %>%
    tibble::enframe(name = NULL)
  colnames(revopvec) <- (c("GenContext.norep"))
  names(prot)[names(prot)=='GenContext.norep'] <- "GenContext.norep.old"
  return(cbind(prot,revopvec))
}

