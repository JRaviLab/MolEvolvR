# Cleanup Genomic Contexts

Cleans up the GenContext column of a data frame by removing certain
characters and rows.

This function removes empty rows based on the 'GenContext' column. A
cleaned up version of the data table is returned.

## Usage

``` r
cleanGenomicContext(
  prot,
  domains_rename = data.frame(old = character(0), new = character(0), stringsAsFactors =
    F),
  condenseRepeatedDomains = TRUE,
  remove_asterisk = TRUE
)
```

## Arguments

- prot:

  A data frame that contains columns 'GenContext.orig'

- domains_rename:

  A data frame containing the domain names to be replaced in a column
  'old' and the replacement in a column 'new'. Defaults to an empty data
  frame with a new and old column such that non of the domains will be
  renamed

- condenseRepeatedDomains:

  Boolean. If TRUE, repeated domains in 'GenContext' are condensed.
  Default is TRUE.

- remove_asterisk:

  Boolean. If TRUE, asterisks in 'ClustName' are removed. Default is
  TRUE.

## Value

A cleaned up version of the data table is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
cleanGenomicContext(prot, domains_rename, T, F)
} # }
```
