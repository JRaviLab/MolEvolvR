# Using the table returned from createIPRScanDomainTable, construct a domain fasta for a single accession number in the original fasta (i.e., the original fasta argument to createIPRScanDomainTable())

Using the table returned from createIPRScanDomainTable, construct a
domain fasta for a single accession number in the original fasta (i.e.,
the original fasta argument to createIPRScanDomainTable())

## Usage

``` r
convertIPRScanDomainTable2FA(df_iprscan_domains)
```

## Arguments

- df_iprscan_domains:

  tbl_df return value from createIPRScanDomainTable

## Value

[Biostrings::AAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
A domain fasta containing all the domains for a single protein in the
original fasta passed as an argument to createIPRScanDomainTable()

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
fasta_domains <- df_iprscan_domains |> convertIPRScanDomainTable2FA()
} # }
```
