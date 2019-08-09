
# operon straightner function #-------

reveql=function(x){
  w=x

  y=rep(NA, length(w))

  d=1

  b=grep("\\*",w)

  for(j in b:length(w)){
    if(w[j]=="=") d= d*(-1)

    if(d==1 && w[j] != "=") {y[j]=paste( w[j], "->",sep = "")

    } else if (d==-1 && w[j] != "="){
      y[j]=paste("<-",w[j], sep = "")

    } else{
      y[j]="="

    }

  }

  if(b >1){
    d=1

  for(j in (b-1):1){
    if(w[j]=="=") d= d*(-1)

    if(d==1 && w[j] != "=") {y[j]=paste( w[j], "->",sep = "")

    } else if (d==-1 && w[j] != "="){
      y[j]=paste("<-",w[j], sep = "")

    } else{
      y[j]="="

    }

  }

  }

  return(y)

}

#----------

 

reverter=function(x){
  opvec=x

  opvec=gsub(pattern = ">",replacement = ">|",x = opvec)

  opvec=gsub(pattern = "<",replacement = "|<",x = opvec)

  opvec=gsub(pattern = "\\|\\|",replacement = "\\|=\\|",x = opvec)

 

  op.list=strsplit(x = opvec, split =  "\\|")

  if(any(is.na(op.list))) op.list[[which(is.na(op.list))]]="-"

  op.list=lapply(1:length(op.list), function(x) {if(any(op.list[[x]]=="")) op.list[[x]][which(op.list[[x]] !="")] else op.list[[x]] })

 

  te=lapply(1:length(op.list), function(x) op.list[[x]][grep("\\*", op.list[[x]])])

  ye=unlist(lapply(te, function(x) substr(x[1],1,1)))

  torev=which(ye=="<")

 

  te=op.list[torev]

  te=lapply(te, function(x) gsub(pattern = "<-|->", replacement = "",x = x))

  te=lapply(te, rev)

  witheq=grep(pattern = "=",x = te)

  withouteq=which(!((1:length(te)) %in% witheq))

  ge=te[witheq]

 

  ge=lapply(1:length(ge), function(x) reveql(ge[[x]]))

  ye=te[withouteq]

  ye=lapply(1:length(ye), function(x) unname(sapply(ye[[x]], function(y) paste(y,"->", sep = ""))))

 

  te[witheq]=ge

  te[withouteq]=ye

  op.list[torev]=te

 

  revopvec=unlist(lapply(op.list, function(x) paste(x, collapse = "")))

  revopvec=gsub(pattern = "=",replacement = "\\|\\|",revopvec)

  return(revopvec)

}

 

##############

colnames(ternary)=c("acc","operon","len", "gen.name","tax","species")

#straighten operons

ternary$operon=reverter(ternary$operon)
