# Lineage Domain Repeats Plot

Lineage Domain Repeats Plot

## Usage

``` r
plotLineageDomainRepeats(query_data, colname, query_DAdoms)
```

## Arguments

- query_data:

  Data frame containing protein homolog data, including relevant domain
  architectures and lineages.

- colname:

  Character. The name of the column in query_data that contains domain
  architectures or other structural information.

- query_DAdoms:

  List or data frame containing 'domains' element, character vector of
  domain names. (NEW ARGUMENT)

## Value

A ggplot object representing a heatmap (tile plot) of domain repeat
counts across different lineages, with color intensity representing the
occurrence of domains.

## Examples

``` r
if (FALSE) { # \dontrun{
plotLineageDomainRepeats()
} # }
```
