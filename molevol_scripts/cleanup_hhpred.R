# Script to analyze hhpred output files

## Dependencies
library(tidyverse); library(here); library(rlang)
library(data.table)
library(d3heatmap) # https://github.com/rstudio/d3heatmap
library(heatmaply) # https://github.com/talgalili/heatmaply

source(here("molevol_scripts/colnames_molevol.R"))
source(here("molevol_scripts/combine_files.R"))

## FILEPATHS
## Assuming that we are starting with "molevol_scripts.Rproj"
gen_path <- here("../molevol_data/common_data/genomes/")
inpath <- here("../molevol_data/project_data/slps/domarch/")

hh_combnd <- combine_files(inpath=inpath, pattern="*.hhr", delim="  ",
                           skip=8, col_names=T)
