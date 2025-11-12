# For a given accession number, get the domain sequences using a interproscan output table & the original FASTA file

For a given accession number, get the domain sequences using a
interproscan output table & the original FASTA file

## Usage

``` r
createIPRScanDomainTable(
  accnum,
  fasta,
  df_iprscan,
  analysis = c("Pfam", "Gene3D")
)
```

## Arguments

- accnum:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  a *single* accession number from the original fasta (fasta param)
  which will be used to search for its sequence's domains (df_iprscan
  param)

- fasta:

  [Biostrings::AAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  original fasta file which was fed into interproscan

- df_iprscan:

  tbl_df the output TSV of interproscan, read as a tibble with
  readIPRScanTSV()

- analysis:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  the domain databases to extract sequences from

## Value

tbl_df table with each domain sequence and a new identifier column

## Examples

``` r
if (FALSE) { # \dontrun{
path_molevol_scripts <- file.path(Sys.getenv("DEV", unset = "/data/molevolvr_transfer/molevolvr_dev"), "molevol_scripts")
setwd(path_molevol_scripts)
source("R/fa2domain.R")
fasta <- Biostrings::readAAStringSet("./tests/example_protein.fa")
df_iprscan <- readIPRScanTSV("./tests/example_iprscan_valid.tsv")
accnum <- df_iprscan$AccNum[1]
df_iprscan_domains <- createIPRScanDomainTable(accnum, fasta, df_iprscan)
} # }
```
