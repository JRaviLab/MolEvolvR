# Cleanup DomArch

Cleanup Domain Architectures

Cleans the DomArch column by replacing/removing certain domains

This function cleans the DomArch column of one data frame by renaming
certain domains according to a second data frame. Certain domains can be
removed according to an additional data frame. The original data frame
is returned with the clean DomArchs column and the old domains in the
DomArchs.old column.

## Usage

``` r
cleanDomainArchitecture(
  prot,
  old = "DomArch.orig",
  new = "DomArch",
  domains_keep,
  domains_rename,
  condenseRepeatedDomains = TRUE,
  removeTails = FALSE,
  removeEmptyRows = F,
  domains_ignore = NULL
)
```

## Arguments

- prot:

  A data frame containing a 'DomArch' column

- old:

  The name of the original column containing domain architecture.
  Defaults to "DomArch.orig".

- new:

  The name of the cleaned column to be created. Defaults to "DomArch".

- domains_keep:

  A data frame containing the domain names to be retained.

- domains_rename:

  A data frame containing the domain names to be replaced in a column
  'old' and the corresponding replacement values in a column 'new'.

- condenseRepeatedDomains:

  Boolean. If TRUE, repeated domains in 'DomArch' are condensed. Default
  is TRUE.

- removeTails:

  Boolean. If TRUE, 'ClustName' will be filtered based on domains to
  keep/remove. Default is FALSE.

- removeEmptyRows:

  Boolean. If TRUE, rows with empty/unnecessary values in 'DomArch' are
  removed. Default is FALSE.

- domains_ignore:

  A data frame containing the domain names to be removed in a column
  called 'domains'

## Value

The original data frame is returned with the clean DomArchs column and
the old domains in the DomArchs.old column.

## Examples

``` r
if (FALSE) { # \dontrun{
cleanDomainArchitecture(prot, TRUE, FALSE,
omains_keep, domains_rename, domains_ignore = NULL)
} # }
```
