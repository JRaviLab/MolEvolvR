# createLineageLookup

Create a look up table that goes from TaxID, to Lineage

## Usage

``` r
createLineageLookup(
  lineage_file = here("data/rankedlineage.dmp"),
  outfile,
  taxonomic_rank = "phylum"
)
```

## Arguments

- lineage_file:

  Path to the rankedlineage.dmp file containing taxid's and their
  corresponding taxonomic rank. rankedlineage.dmp can be downloaded at
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

- outfile:

  File the resulting lineage lookup table should be written to

- taxonomic_rank:

  The upperbound of taxonomic rank that the lineage includes. The
  lineaege will include superkingdom\>...\>taxonomic_rank. Choices
  include: "supperkingdom", "phylum", "class","order", "family",
  "genus", and "species"

## Value

A tibble containing the tax IDs and their respective lineages up to the
specified taxonomic rank, saved as a tab-separated file.

## Author

Samuel Chen

## Examples

``` r
if (FALSE) { # \dontrun{
createLineageLookup(lineage_file = "data/rankedlineage.dmp",
                     outfile = "data/lineage_lookup.tsv",
                     taxonomic_rank = "family")
} # }
```
