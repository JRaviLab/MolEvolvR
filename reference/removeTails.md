# Remove Tails

Remove tails/singletons

This function ... Certain low frequency domain architectures can be
removed. The original data frame is returned with the corresponding
cleaned up column.

## Usage

``` r
removeTails(prot, by_column = "DomArch", keep_domains = FALSE)
```

## Arguments

- prot:

  A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.

- by_column:

  Default column is 'DomArch'. Can also take 'ClustName', 'GenContext'
  as input.

- keep_domains:

  Default is False Keeps tail entries that contain the query domains.

## Value

The original data frame with singletons removed from the specified
column.

## Examples

``` r
if (FALSE) { # \dontrun{
removeTails(prot, "DomArch")
} # }
```
