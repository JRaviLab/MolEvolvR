##reading domains.rename
#Can just use psp again
#test run on SIg+TM+LADB... or DomArch
library(tidyverse)
domains_rename <- read_delim("C:/Users/samue/Desktop/domains.rename", "\t",

                             escape_double = FALSE, col_names = FALSE,

                             trim_ws = TRUE)
#data must have 2 columns
names(domains_rename)=c("old","new") #rename one column to old and the other to new

domains_rename

eDNA_data$archold=eDNA_data$arch
#eDNA is cleanblastclust file, probably could test on DomArch

#For loop iterates through the rows from 1 to the end of domains_rename$old
for (j in 1:length(domains_rename$old)) {

  #eDNA_data$arch[grep(domains_rename$old[j], eDNA_data$arch)]=gsub(domains_rename$old[j], domains_rename$new[j],x = eDNA_data$arch[grep(domains_rename$old[j], eDNA_data$arch)])
#regterm is the string "^domainofold\\+"
  regterm=paste ("^",domains_rename$old[j],"\\+",sep="")
#repterm is the string "domofnew+"
  repterm=paste (domains_rename$new[j],"\\+",sep="")
#All indeces in eDNA_data$arch that match regterm, are replaced with repterm
  eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

#regterm now is "+domainofold$" (very similar to first regterm)
  regterm=paste ("\\+",domains_rename$old[j],"$",sep="")
#repterm is now "+domainofnew" (very similar to first repterm)
  repterm=paste ("\\+",domains_rename$new[j],sep="")
#All indeces in eDNA_data$arch that match second regterm are replaced with second repterm
  eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

#regterm3 = "+domainofold+"
  regterm=paste ("\\+",domains_rename$old[j],"\\+",sep="")
#repterm3 = "+domainofnew+"
  repterm=paste ("\\+",domains_rename$new[j],"\\+",sep="")

  #repterm= domains_rename$new[j]
#What does[1] do?
#while iterates keeps going until no more strings in eDNA_data$arch match regterm
  while ( length(eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]) > 0) {
#replace the regterm with repterm
    eDNA_data$arch[grep(regterm, eDNA_data$arch)]=gsub(regterm, repterm,x = eDNA_data$arch[grep(regterm, eDNA_data$arch)])

  }

}



##reading domains.afterrename.group.tsv

domains_afterrename_group <- read_delim("C:/Users/samue/Desktop/domains.afterrename.group.tsv",

                                        "\t", escape_double = FALSE, col_names = FALSE,

                                        trim_ws = TRUE)

domains_afterrename_group





#generating domains to ignore and cleaning the names in arch
#ge is a vector containing vectors of the strings in eDNA_data$arch that are split where there is a '\+'
ge = unlist(strsplit(eDNA_data$arch, "\\+"))

#domains.ignore = the elements in ge that aren't in domains_afterrename_group$X1
domains.ignore=setdiff(ge, domains_afterrename_group$X1)
#for loop essentially replaces the elements in eDNA_data$arch that match domains.ignore with
#\+, or maybe just delete them all
#for loop goes from 1 to length of domains.ignore
for (j in 1:length(domains.ignore)) {
#regterm = "\+domains.ignore\+", changes for each iteration
  regterm=paste ("\\+",domains.ignore[j],"\\+",sep="")
#while loop goes until no more terms in eDNA_data$arch match regterm
  while ( length(eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]) > 0) {
#substitute "\+" for regterm in eDNA_data$arch
   eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]=gsub(regterm[1], "\\+",x = eDNA_data$arch[grep(regterm[1], eDNA_data$arch)] )

  }
#regterm = "^domains.ignore\+|\+domains.ignore$"
  regterm=paste ("^",domains.ignore[j],"\\+","|","\\+",domains.ignore[j],"$",sep="")
#replace eDNA_data$arch elements that match regterm with "" (delete the elements that match domains.ignore)
  eDNA_data$arch[grep(regterm[1], eDNA_data$arch)]=gsub(regterm[1], "",x = eDNA_data$arch[grep(regterm[1], eDNA_data$arch)] )



  #eDNA_data$arch=gsub("\\+{2,}", "+", eDNA_data$arch)

  #eDNA_data$arch=gsub("^\\+|\\+$", "", eDNA_data$arch)

}



#cleaning up orphan SIG/TM archs

eDNA_data=eDNA_data[- which (eDNA_data$arch == "SIG" | eDNA_data$arch == "TM"), ]
