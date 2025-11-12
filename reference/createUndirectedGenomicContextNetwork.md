# createUndirectedGenomicContextNetwork

This function creates a domain network from the 'DomArch' column.

A network of domains is returned based on shared domain architectures.

## Usage

``` r
createUndirectedGenomicContextNetwork(
  prot,
  column = "GenContext",
  domains_of_interest,
  cutoff_type = "Lineage",
  cutoff = 1,
  layout = "grid"
)
```

## Arguments

- prot:

  A data frame that contains the column 'DomArch'.

- column:

  Name of column containing Domain architecture from which nodes and
  edges are generated.

- domains_of_interest:

  Character vector specifying the domains of interest.

- cutoff_type:

  Character. Used to determine how data should be filtered. Either

  - "Total Count" to filter off the total amount of times a domain
    architecture occurs

- cutoff:

  Integer. Only use domains that occur at or above the cutoff for total
  counts if cutoff_type is "Total Count". Only use domains that appear
  in cutoff or greater lineages if cutoff_type is Lineage.

- layout:

  Character. Layout type to be used for the network. Options are:

  - "circle"

  - "random"

  - "auto"

## Value

A plot of the domain architecture network.

## Examples

``` r
if (FALSE) { # \dontrun{
createUndirectedGenomicContextNetwork(pspa, column = "DomArch",
domains_of_interest = c("Domain1", "Domain2"),
cutoff_type = "Total Count", cutoff = 10)
} # }
```
