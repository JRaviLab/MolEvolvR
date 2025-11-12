# Cleanup FASTA Header

Parse accesion numbers from fasta and add a suffix of the ith occurence
to handle duplicates

## Usage

``` r
cleanFAHeaders(fasta)
```

## Arguments

- fasta:

  An
  [Biostrings::XStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object representing the sequences from a FASTA file. The sequence
  names (headers) will be adjusted for uniqueness and sanitized.

## Value

[Biostrings::XStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
fasta with adjusted names (headers)

## Examples

``` r
if (FALSE) { # \dontrun{
AAStringSet(c("xxx" = "ATCG", "xxx" = "GGGC")) |>
    cleanFAHeaders()
} # }
```
