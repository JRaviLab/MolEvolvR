# Remove Empty

Remove empty rows by column

Removes empty rows in the specified column.

This function ... The original data frame is returned with the
corresponding cleaned up column.

## Usage

``` r
removeEmptyRows(prot, by_column = "DomArch")
```

## Arguments

- prot:

  A data frame containing 'DomArch', 'Species', 'GenContext',
  'ClustName' columns.

- by_column:

  Column by which empty rows should be removed to domain+domain -\>
  domain(s). Default column is 'DomArch'. Can also take the following as
  input, 'Species', 'GenContext', 'ClustName'.

## Value

A tibble with rows removed where the specified column contains `"-"`,
`"NA"`, or an empty string.

## Examples

``` r
if (FALSE) { # \dontrun{
removeEmptyRows(prot, "DomArch")
} # }
```
