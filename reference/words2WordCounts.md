# words2WordCounts

Get word counts (wc) DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)

## Usage

``` r
words2WordCounts(string)
```

## Arguments

- string:

  A character string containing the elements (words) to count. This
  would typically be a space-delimited string representing domain
  architectures or genomic contexts.

## Value

A tibble (tbl_df) with two columns:

- `words`:

  A column containing the individual words (domains or domain
  architectures).

- `freq`:

  A column containing the frequency counts for each word.

## Examples

``` r
if (FALSE) { # \dontrun{
tibble::tibble(DomArch = c("aaa+bbb", "a+b", "b+c", "b-c")) |>
    elements2Words() |>
    words2WordCounts()
} # }
```
