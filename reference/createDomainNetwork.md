# Domain Network

This function creates a domain network from the 'DomArch' column.

A network of domains is returned based on shared domain architectures.

## Usage

``` r
createDomainNetwork(
  prot,
  column = "DomArch",
  domains_of_interest,
  cutoff = 70,
  layout = "nice",
  query_color = adjustcolor("green", alpha.f = 0.5)
)
```

## Arguments

- prot:

  A data frame that contains the column 'DomArch'.

- column:

  Name of column containing Domain architecture from which nodes and
  edges are generated.

- domains_of_interest:

  Character vector specifying domains of interest.

- cutoff:

  Integer. Only use domains that occur at or above the cutoff for total
  counts if cutoff_type is "Total Count". Only use domains that appear
  in cutoff or greater lineages if cutoff_type is Lineage.

- layout:

  Character. Layout type to be used for the network. Options are:

  - "circle"

  - "random"

  - "auto"

- query_color:

  Character. Color to represent the queried domain in the network.

## Value

A network visualization of domain architectures.

## Examples

``` r
if (FALSE) { # \dontrun{
createDomainNetwork(pspa)
} # }
```
