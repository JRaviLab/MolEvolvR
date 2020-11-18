# SCRIPT to Cleanup VFDB file and convert it to 'Tidy' format
# Created: Mar 14, 2019 | Janani Ravi

##################
# library(tidyverse)
# library(here)
##################

## Suggested Reading/Reference material
## 1. Books: R4DS, HOPR
## 2. Tidyverse workshop material from RLEL, Feb 2019
## 3. ? or ??FUNCTION_NAME within R to check what the function does
## 4. Cheatsheets for dplyr, readr, basic regex (regular expression)

## Data Import
here() # here's where R is going to look for your data file
# Reading the raw text file from VFDB --> Mycobacterium
mycobacteria_vfs <- read_tsv("data/Mycobacterium_VFs_comparsion.txt",
                             skip=1) # skipping the first row (not the column header row)

## LOOKING at your data
# Checking to see if it looks alright --> column names look too complicated!
# View(mycobacteria_vfs) # opens in a separate window
glimpse(mycobacteria_vfs)

## Changing column names
# first storing the current column names in a dataframe/tibble with the column name, "old"
colnames_myco_vfs <- tibble(old=colnames(mycobacteria_vfs))

## Mission: Cleanup!! Using mutate & gsub
# Using mutate from tidyverse/dplyr; Creating a new variable "new"
# Replacing space, dot, brackets and forward/backslash with underscore "_"
# Note: Check after each mutate step to see what's happening
colnames_myco_vfs <- colnames_myco_vfs %>%
  mutate(new = gsub(pattern="[ \\.\\(\\)-///]", # // to release the special character
                    replacement="_",
                    x=old)) %>%
  # then replacing multiple underscores w/ a single underscore
  mutate(new = gsub(pattern="[\\_]+", # plus is for repeated characters
                    replacement="_",
                    x=new)) %>%
  # to remove 'chromosome_' because that makes no sense in this Mycobacterial context
  mutate(new = gsub("chromosome_", "", new)) %>%
  # removing length at the end XYZ_bp
  mutate(new = gsub("_[0-9]+_bp_$",
                    "",
                    new)) %>%
  # removing extra _ after NC
  mutate(new = gsub("_NC_",
                    "_NC",
                    new))

## Now replace the column names of the original VFDB downloaded txt file w/ the cleaned up version
colnames(mycobacteria_vfs) <- colnames_myco_vfs$new
# Check again to see if everything looks alright
glimpse(mycobacteria_vfs)


## Convert WIDE -> LONG (tidy data) format?
mycobacteria_vfs_long <- mycobacteria_vfs %>%
  gather(-VFclass, -Virulence_factors, -Related_genes,
         matches("^M_"), key=species, value=gene_name)

# Check again!
glimpse(mycobacteria_vfs_long)

## Saving your data
here() # here's where your data is going to be stored!
write_tsv(x=mycobacteria_vfs,
          path="data/mycobacteria_vfs_wide.txt", col_names=T)
write_tsv(x=mycobacteria_vfs_long,
          path="data/mycobacteria_vfs_long.txt", col_names=T)

## On to other comparative pathogenomics analyses to identify genes unique to all avium species or all pathogenic NTM or all M. tb complex genomes etc ...
