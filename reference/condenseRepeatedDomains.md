# condenseRepeatedDomains

Condense repeated domains

Condenses repeated domains in the specified column.

This function identifies repeated domains and condenses them to (s). ??
Certain domains can be removed according to an additional data frame.
The original data frame is returned with the corresponding cleaned up
column.

## Usage

``` r
condenseRepeatedDomains(prot, by_column = "DomArch", excluded_prots = c())
```

## Arguments

- prot:

  A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.

- by_column:

  Column in which repeats are condensed to domain+domain -\> domain(s).

- excluded_prots:

  Vector of strings that condenseRepeatedDomains should not reduce to
  (s). Defaults to c()

## Value

A data frame with condensed repeated domains in the specified column.

## Examples

``` r
if (FALSE) { # \dontrun{
condenseRepeatedDomains(prot, "DomArch")
} # }
```
