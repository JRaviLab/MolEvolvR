# prepareColumnParams

prepareColumnParams

## Usage

``` r
prepareColumnParams(count_data, fill_by_n, sort_by_n)
```

## Arguments

- count_data:

  A data frame containing the data.

- fill_by_n:

  Logical indicating if fill color is based on counts.

- sort_by_n:

  Logical indicating if data should be sorted by counts.

## Value

A data frame of parameters for treemap visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
count_data <- data.frame(Category = c("A", "B", "C"),
                          n = c(10, 20, 15))
params <- prepareColumnParams(count_data, fill_by_n = TRUE, sort_by_n = FALSE)
params
} # }
```
