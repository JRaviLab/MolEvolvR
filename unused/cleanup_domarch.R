## Functions to cleanup the domain architecture column
## Originally from V&L
## Modified by JR, SZC.

#################
## Pkgs needed ##
#################
suppressPackageStartupMessages(library(tidyverse))

## Reading domains.rename
# Can just use psp again
# Test run on SIg+TM+LADB... or DomArch
domains_rename <- read_tsv("data/acc_files/domains_rename.txt")
# data must have 2 columns (old, new)
# names(domains_rename) <- c("old", "new") #rename one column to old and the other to new

domains_rename

## Reading input protein/domain file
prot <- read_tsv("data/rawdata_tsv/all_raw.txt")

# edited to create a duplicate 'new' DomArch column
prot$DomArch <- prot$DomArch.orig
#prot/eDNA_data is cleanblastclust file, probably could test on DomArch

# For loop iterates through the rows from 1 to the end of domains_rename$old
for (j in 1:length(domains_rename$old)) {
  # ?? Replace globally? Or comment out??
  # prot$DomArch.orig[grep(pattern=domains_rename$old[j],
  #                x=prot$DomArch.orig)] <- gsub(pattern=domains_rename$old[j],
  #                                         replacement=domains_rename$new[j],
  #                                         x=prot$DomArch.orig[grep(domains_rename$old[j],
  #                                                                  prot$DomArch.orig)])
  # Tidyverse version
  # prot$DomArch <- prot$DomArch.orig %>%
  #   str_replace_all(pattern=domains_rename$old[j], replacement=domains_rename$new[j])

  ## Replace first instance
  #find_term #1: "^domain_old\\+"
  find_term <- paste("^", domains_rename$old[j],
                   "\\+", sep="")
  #replace_term #1: "domain_new+"
  replace_term <- paste(domains_rename$new[j],
                   "\\+", sep="")
  #All indices in prot$DomArch that match find_term, are replaced with replace_term
  prot$DomArch[grep(find_term, prot$DomArch)] <- gsub(find_term, replace_term,
                                                      x=prot$DomArch[grep(find_term,
                                                                          prot$DomArch)])

  ## Replace last instance
  #find_term #2: "+domain_old$" (very similar to first find_term)
  find_term <- paste("\\+", domains_rename$old[j],
                   "$", sep="")
  #replace_term #2: "+domain_new" (very similar to first replace_term)
  replace_term <- paste("\\+", domains_rename$new[j],
                   sep="")
  #All indices in prot$DomArch that match second find_term are replaced with second replace_term
  prot$DomArch[grep(find_term, prot$DomArch)] <- gsub(find_term, replace_term,
                                                      x=prot$DomArch[grep(find_term,
                                                                          prot$DomArch)])

  ## Replace middle instances
  #find_term #3: "+domain_old+"
  find_term <- paste("\\+", domains_rename$old[j],
                   "\\+", sep="")
  #replace_term #3: "+domainofnew+"
  replace_term <- paste("\\+", domains_rename$new[j],
                   "\\+", sep="")

  #replace_term= domains_rename$new[j]
  #What does[1] do?
  #while iterates keeps going until no more strings in prot$DomArch match find_term
  while(length(prot$DomArch[grep(find_term[1], prot$DomArch)]) > 0) {
    #replace the find_term with replace_term
    prot$DomArch[grep(find_term, prot$DomArch)] <- gsub(find_term, replace_term,
                                                        x=prot$DomArch[grep(find_term,
                                                                            prot$DomArch)])
  }
}

## UNUSED by JR & co. so far!
## Generating ignore list! Reading in our list instead.
# ## Reading domains.afterrename.group.tsv
#
# domains_afterrename_group <- read_delim("data/acc_files/domains.afterrename.group.tsv",
#                                         "\t", escape_double = FALSE, col_names = FALSE,
#                                         trim_ws = TRUE)
# domains_afterrename_group
#
#
# # Generating domains to ignore and cleaning the names in arch
# #ge is a vector containing vectors of the strings in prot$DomArch that are split where there is a '\+'
# ge <- unlist(strsplit(prot$DomArch, "\\+"))
#
# #domains_ignore = the elements in ge that aren't in domains_afterrename_group$X1
# domains_ignore <- setdiff(ge, domains_afterrename_group$X1)


## Reading domains_ignore
domains_ignore <- read_tsv("data/acc_files/domains_ignore.txt")
# for loop essentially replaces the elements in prot$DomArch that match domains_ignore with
#\+, or maybe just delete them all

# for loop goes from 1 to length of domains_ignore
for (j in 1:length(domains_ignore)) {
  #find_term = "\+domains_ignore\+", changes for each iteration
  find_term <- paste("\\+", domains_ignore[j], "\\+",
                   sep="")
  #while loop goes until no more terms in prot$DomArch match find_term
  while(length(prot$DomArch[grep(find_term[1], prot$DomArch)]) > 0) {
    #substitute "\+" for find_term in prot$DomArch
    prot$DomArch[grep(find_term[1], prot$DomArch)] <- gsub(find_term[1], "\\+",
                                                           x=prot$DomArch[grep(find_term[1],
                                                                               prot$DomArch)])

  }
  #find_term = "^domains_ignore\+|\+domains_ignore$"
  find_term <- paste("^", domains_ignore[j], "\\+", "|", "\\+", domains_ignore[j], "$",
                   sep="")
  #replace prot$DomArch elements that match find_term with "" (delete the elements that match domains_ignore)
  prot$DomArch[grep(find_term[1], prot$DomArch)] <- gsub(find_term[1], "",
                                                         x=prot$DomArch[grep(find_term[1],
                                                                             prot$DomArch)] )

  #prot$DomArch=gsub("\\+{2, }", "+", prot$DomArch)
  #prot$DomArch=gsub("^\\+|\\+$", "", prot$DomArch)
}


## Cleaning up orphan SIG/TM archs
prot <- prot[- which (prot$DomArch == "SIG" | prot$DomArch == "TM"), ]
