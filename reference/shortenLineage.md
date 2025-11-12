# shortenLineage

This function abbreviates lineage names by shortening the first part of
the string (up to a given delimiter).

## Usage

``` r
shortenLineage(data, colname = "Lineage", abr_len = 1)
```

## Arguments

- data:

  A data frame that contains a column with lineage names to be
  shortened.

- colname:

  Character. The name of the column in the data frame containing the
  lineage strings to be shortened. Default is `"Lineage"`.

- abr_len:

  Integer. The number of characters to retain after the first letter. If
  set to 1, only the first letter of each segment before the delimiter
  (`>`) is retained. Default is 1.

## Value

A modified data frame where the specified lineage column has been
shortened.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(Lineage = c("Bacteria>Firmicutes>Clostridia",
"Archaea>Euryarchaeota>Thermococci"))
shortened_df <- shortenLineage(df, colname = "Lineage", abr_len = 1)
shortened_df
} # }
```
