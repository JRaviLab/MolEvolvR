# Cleanup Species

Cleanup Species

Cleans up the species column of a data frame by removing certain
characters and rows.

This function removes unneccessary characters from the 'Species' column.
A cleaned up version of the data table is returned.

## Usage

``` r
cleanSpecies(prot, removeEmptyRows = FALSE)
```

## Arguments

- prot:

  A data frame that contains columns 'Species'.

- removeEmptyRows:

  Boolean. If TRUE, rows with empty/unnecessary values in 'Species' are
  removed. Default is false.

## Value

The original data frame with Species cleaned.

## Examples

``` r
if (FALSE) { # \dontrun{
cleanSpecies(prot, TRUE)
} # }
```
