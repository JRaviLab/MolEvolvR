# countByColumn

Function to obtain element counts (DA, GC)

## Usage

``` r
countByColumn(prot = prot, column = "DomArch", min.freq = 1)
```

## Arguments

- prot:

  A data frame containing the dataset to analyze, typically with
  multiple columns including the one specified by the `column`
  parameter.

- column:

  A character string specifying the name of the column to analyze. The
  default is "DomArch".

- min.freq:

  An integer specifying the minimum frequency an element must have to be
  included in the output. Default is 1.

## Value

A tibble with two columns:

- `column`:

  The unique elements from the specified column (e.g., "DomArch").

- `freq`:

  The frequency of each element, i.e., the number of times each element
  appears in the specified column.

The tibble is filtered to only include elements that have a frequency
greater than or equal to `min.freq` and does not include elements with
`NA` values or those starting with a hyphen ("-").

## Examples

``` r
if (FALSE) { # \dontrun{
countByColumn(prot = my_data, column = "DomArch", min.freq = 10)
} # }
```
