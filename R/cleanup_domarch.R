##reading domains.rename

domains_rename <- read_delim("domains.rename", "\t",

                             escape_double = FALSE, col_names = FALSE,

                             trim_ws = TRUE)

names(domains_rename)=c("old","new")

domains_rename

eDNA_data$archold=eDNA_data$arch

for (j in 1:length(domains_rename$old)) {

  #eDNA_data$arch[grep(domains_rename$old[j], eDNA_data$arch)]=gsub(domains_rename$old[j], domains_rename$new[j],x = eDNA_data$arch[grep(domains_rename$old[j], eDNA_data$arch)])

  regterm=paste ("^",domains_rename$old[j],"\\+",sep="")

  repterm=paste (domains_rename$new[j],"\\+",sep="")

  eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

 

  regterm=paste ("\\+",domains_rename$old[j],"$",sep="")

  repterm=paste ("\\+",domains_rename$new[j],sep="")

  eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

 

  regterm=paste ("\\+",domains_rename$old[j],"\\+",sep="")

  repterm=paste ("\\+",domains_rename$new[j],"\\+",sep="")

  #repterm= domains_rename$new[j]

  while ( length(eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]) > 0) {

    eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

  }

}

 

##reading domains.afterrename.group.tsv

domains_afterrename_group <- read_delim("domains.afterrename.group.tsv",

                                        "\t", escape_double = FALSE, col_names = FALSE,

                                        trim_ws = TRUE)

domains_afterrename_group

 

 

#generating domains to ignore and cleaning the names in arch

ge = unlist(strsplit(eDNA_data$arch, "\\+"))

domains.ignore=setdiff(ge, domains_afterrename_group$X1)

for (j in 1:length(domains.ignore)) {

  regterm=paste ("\\+",domains.ignore[j],"\\+",sep="")

  while ( length(eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]) > 0) {

   eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]=gsub(regterm[1], "\\+",x = eDNA_data$arch[grep(regterm[1], eDNA_data$arch)] )

  }

  regterm=paste ("^",domains.ignore[j],"\\+","|","\\+",domains.ignore[j],"$",sep="")

  eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]=gsub(regterm[1], "",x = eDNA_data$arch[grep(regterm[1], eDNA_data$arch)] )

 

  #eDNA_data$arch=gsub("\\+{2,}", "+", eDNA_data$arch)

  #eDNA_data$arch=gsub("^\\+|\\+$", "", eDNA_data$arch)

}

 

#cleaning up orphan SIG/TM archs

eDNA_data=eDNA_data[- which (eDNA_data$arch == "SIG" | eDNA_data$arch == "TM"), ]
