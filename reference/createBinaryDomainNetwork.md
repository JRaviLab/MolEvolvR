# Domain Network

This function creates a domain network from the 'DomArch' column.

Domains that are part of the 'domains_of_interest' are a different node
color than the other domains.

A network of domains is returned based on shared domain architectures.

## Usage

``` r
createBinaryDomainNetwork(
  prot,
  column = "DomArch",
  domains_of_interest,
  cutoff = 70,
  layout = "nice",
  query_color = adjustcolor("yellow", alpha.f = 0.5),
  partner_color = adjustcolor("skyblue", alpha.f = 0.5),
  border_color = adjustcolor("grey", alpha.f = 0.8),
  IsDirected = T
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

  Color that the nodes of the domains in the domains_of_interest vector
  are colored

- partner_color:

  Color that the nodes that are not part of the domains_of_interest
  vector are colored

- border_color:

  Color for the borders of the nodes.

- IsDirected:

  Is the network directed? Set to false to eliminate arrows

## Value

A network visualization of domain architectures.

## Examples

``` r
if (FALSE) { # \dontrun{
createBinaryDomainNetwork(pspa)
} # }
```
