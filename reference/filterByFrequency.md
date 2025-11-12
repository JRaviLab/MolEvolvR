# filterByFrequency

Function to filter based on frequencies

## Usage

``` r
filterByFrequency(x, min.freq)
```

## Arguments

- x:

  A tibble (tbl_df) containing at least two columns: one for elements
  (e.g., `words`) and one for their frequency (e.g., `freq`).

- min.freq:

  A numeric value specifying the minimum frequency threshold. Only
  elements with frequencies greater than or equal to this value will be
  retained.

## Value

A tibble with the same structure as `x`, but filtered to include only
rows where the frequency is greater than or equal to `min.freq`.

## Examples

``` r
if (FALSE) { # \dontrun{
filterByFrequency()
} # }
```
