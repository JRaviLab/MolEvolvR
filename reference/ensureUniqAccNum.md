# make accnums unique

Append an index of occurence suffix to each accession number (or any
character vector) making them unique

## Usage

``` r
ensureUniqAccNum(accnums)
```

## Arguments

- accnums:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  a vector of accession numbers

## Value

[rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
a vector of adjusted, unique accession numbers

## Examples

``` r
if (FALSE) { # \dontrun{
c("xxx", "xxx", "xxx", "yyy", "yyy") |>
    ensureUniqAccNum()
} # }
```
