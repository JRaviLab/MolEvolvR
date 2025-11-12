# getDomainsFromFA

getDomainsFromFA

## Usage

``` r
getDomainsFromFA(
  fasta,
  df_iprscan,
  analysis = c("Pfam", "Gene3D"),
  verbose = FALSE
)
```

## Arguments

- fasta:

  [Biostrings::AAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  a protein (AA) fasta

- df_iprscan:

  tbl_df the interproscan results from the original fasta

- analysis:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  the domain databases to extract sequences from

## Value

fasta_domains
[Biostrings::AAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
fasta of domains

## Examples

``` r
if (FALSE) { # \dontrun{
path_molevol_scripts <- file.path(Sys.getenv("DEV", unset = "/data/molevolvr_transfer/molevolvr_dev"), "molevol_scripts")
setwd(path_molevol_scripts)
source("R/fa2domain.R")
fasta <- Biostrings::readAAStringSet("./tests/example_protein.fa")
df_iprscan <- readIPRScanTSV("./tests/example_iprscan_valid.tsv")
getDomainsFromFA(fasta, df_iprscan)
} # }
```
