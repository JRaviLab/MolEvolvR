# Genomic Context Directed Network

This function creates a Genomic Context network from the 'GenContext'
column.

A network of Genomic Context is returned.

## Usage

``` r
createGenomicContextNetwork(
  prot,
  domains_of_interest,
  column = "GenContext",
  cutoff = 40,
  layout = "grid",
  directed = TRUE
)
```

## Arguments

- prot:

  A data frame that contains the column 'GenContext'.

- domains_of_interest:

  Character vector of domains of interest.

- column:

  Name of column containing Genomic Context from which nodes and edges
  are generated.

- cutoff:

  Integer. Only use GenContexts that occur at or above the cutoff
  percentage for total count

- layout:

  Character. Layout type to be used for the network. Options are:

  - "circle"

  - "random"

  - "auto"

  - "nice"

- directed:

  Is the network directed?

## Value

A plot of the genomic context network.

## Examples

``` r
if (FALSE) { # \dontrun{
gc_directed_network(pspa, column = "GenContext", cutoff = 55)
} # }
```
