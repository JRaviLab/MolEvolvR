suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidytext))

prot <- read_tsv("data/rawdata_tsv/all_raw.txt")

df <- prot$DomArch.orig

temp <- df %>%
  gsub(pattern="\\+", replacement=" ") %>%
  tibble(doms=.)

temp2 <- df %>%
  tibble::as_tibble() %>%
  dplyr::mutate_if(is.character,
                 stringr::str_replace_all,
                 pattern = "\\+",
                 replacement = " ") %>%
  dplyr::mutate_if(is.character,
                   stringr::str_replace_all,
                   pattern = "-",
                   replacement = "XYZ")
temp2 %>%
  unnest_tokens(ngram, value, token = "regex", pattern="\n")

temp2 %>%
  unnest_tokens(output=bigram, input=value, token="ngrams", n=2, format="text")
  #unnest_tokens(bigram, doms, token = stringr::str_split, pattern = " ")

