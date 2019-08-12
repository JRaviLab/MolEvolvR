## Functions to cleanup the domain architecture column

library(tidyverse)

## Reading domains.rename
# Can just use psp again
# Test run on SIg+TM+LADB... or DomArch
domains_rename <- read_delim("data/acc_files/domains.replace.txt",
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE)
# data must have 2 columns (old, new)
names(domains_rename) <- c("old","new") #rename one column to old and the other to new

domains_rename

## Reading input protein/domain file
prot <- read_delim("data/rawdata_tsv/all.txt")

prot$archold <- prot$DomArch
#prot/eDNA_data is cleanblastclust file, probably could test on DomArch

# For loop iterates through the rows from 1 to the end of domains_rename$old
for (j in 1:length(domains_rename$old)) {

  #prot$arch[grep(domains_rename$old[j], prot$arch)]=gsub(domains_rename$old[j], domains_rename$new[j],x = prot$arch[grep(domains_rename$old[j], prot$arch)])
  #regterm is the string "^domainofold\\+"
  regterm <- paste("^",domains_rename$old[j],
                   "\\+",sep="")
  #repterm is the string "domofnew+"
  repterm <- paste(domains_rename$new[j],
                   "\\+",sep="")
  #All indeces in prot$arch that match regterm, are replaced with repterm
  prot$arch[grep(regterm, prot$arch)] <- gsub(regterm, repterm,
                                              x = prot$arch[grep(regterm, prot$arch)])

  #regterm now is "+domainofold$" (very similar to first regterm)
  regterm <- paste("\\+",domains_rename$old[j],
                   "$",sep="")
  #repterm is now "+domainofnew" (very similar to first repterm)
  repterm <- paste("\\+",domains_rename$new[j],
                   sep="")
  #All indeces in prot$arch that match second regterm are replaced with second repterm
  prot$arch[grep(regterm, prot$arch)] <- gsub(regterm, repterm,
                                              x = prot$arch[grep(regterm, prot$arch)])

  #regterm3 = "+domainofold+"
  regterm <- paste("\\+", domains_rename$old[j],
                   "\\+", sep="")
  #repterm3 = "+domainofnew+"
  repterm <- paste("\\+", domains_rename$new[j],
                   "\\+", sep="")

  #repterm= domains_rename$new[j]
  #What does[1] do?
  #while iterates keeps going until no more strings in prot$arch match regterm
  while(length(prot$arch[grep(regterm[1], prot$arch)]) > 0) {
    #replace the regterm with repterm
    prot$arch[grep(regterm, prot$arch)] <- gsub(regterm, repterm,
                                                x = prot$arch[grep(regterm, prot$arch)])
  }
}

## Reading domains.afterrename.group.tsv

domains_afterrename_group <- read_delim("data/acc_files/domains.afterrename.group.tsv",
                                        "\t", escape_double = FALSE, col_names = FALSE,
                                        trim_ws = TRUE)
domains_afterrename_group


# Generating domains to ignore and cleaning the names in arch
#ge is a vector containing vectors of the strings in prot$arch that are split where there is a '\+'
ge <- unlist(strsplit(prot$arch, "\\+"))

#domains.ignore = the elements in ge that aren't in domains_afterrename_group$X1
domains.ignore <- setdiff(ge, domains_afterrename_group$X1)

# for loop essentially replaces the elements in prot$arch that match domains.ignore with
#\+, or maybe just delete them all

# for loop goes from 1 to length of domains.ignore
for (j in 1:length(domains.ignore)) {
  #regterm = "\+domains.ignore\+", changes for each iteration
  regterm <- paste("\\+", domains.ignore[j], "\\+",
                   sep="")
  #while loop goes until no more terms in prot$arch match regterm
  while(length(prot$arch[grep(regterm[1], prot$arch)]) > 0) {
    #substitute "\+" for regterm in prot$arch
    prot$arch[grep(regterm[1], prot$arch)] <- gsub(regterm[1], "\\+",
                                                   x = prot$arch[grep(regterm[1], prot$arch)])

  }
  #regterm = "^domains.ignore\+|\+domains.ignore$"
  regterm <- paste("^", domains.ignore[j], "\\+","|","\\+", domains.ignore[j],"$",
                   sep="")
  #replace prot$arch elements that match regterm with "" (delete the elements that match domains.ignore)
  prot$arch[grep(regterm[1], prot$arch)] <- gsub(regterm[1], "",
                                                 x = prot$arch[grep(regterm[1], prot$arch)] )

  #prot$arch=gsub("\\+{2,}", "+", prot$arch)
  #prot$arch=gsub("^\\+|\\+$", "", prot$arch)
}


## Cleaning up orphan SIG/TM archs
prot <- prot[- which (prot$arch == "SIG" | prot$arch == "TM"), ]
