# Cleanup Clust

Cleanup cluster file

Cleans a cluster file by removing rows that do not contain the query in
the cluster.

This function removes irrelevant rows which do not contain the query
protein within the ClustName column. The return value is the cleaned up
data frame.

## Usage

``` r
cleanClusters(
  prot,
  domains_rename,
  domains_keep,
  condenseRepeatedDomains = TRUE,
  removeTails = FALSE,
  removeEmptyRows = FALSE
)
```

## Arguments

- prot:

  A data frame that must contain columns Query and ClustName.

- domains_rename:

  A data frame containing the domain names to be replaced in a column
  'old' and the corresponding replacement values in a column 'new'.

- domains_keep:

  A data frame containing the domain names to be retained.

- condenseRepeatedDomains:

  Boolean. If TRUE, repeated domains in 'ClustName' are condensed.
  Default is TRUE.

- removeTails:

  Boolean. If TRUE, 'ClustName' will be filtered based on domains to
  keep/remove. Default is FALSE.

- removeEmptyRows:

  Boolean. If TRUE, rows with empty/unnecessary values in 'ClustName'
  are removed. Default is FALSE.

## Value

Cleaned up data frame

## Examples

``` r
if (FALSE) { # \dontrun{
cleanClusters(prot, TRUE, FALSE, domains_keep, domains_rename)
} # }
```
