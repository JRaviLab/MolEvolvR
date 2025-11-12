# prepareSingleColumnParams

prepareSingleColumnParams

## Usage

``` r
prepareSingleColumnParams(df, col_num, root)
```

## Arguments

- df:

  A data frame containing the data to be processed.

- col_num:

  An integer representing the column number to process.

- root:

  A string representing the root node for the treemap.

## Value

A data frame containing parameters for the specified column for treemap
visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(Category = c("A", "A", "B", "B", "C"),
                 n = c(10, 20, 30, 40, 50))
params <- prepareSingleColumnParams(df, col_num = 1, root = "Root")
params
} # }
```
