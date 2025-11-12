# Combining clean ipr files

Combining clean ipr files

## Usage

``` r
combineIPR(inpath, ret = FALSE)
```

## Arguments

- inpath:

  Character. The path to the directory containing the `.iprscan_cln.tsv`
  files to be combined.

- ret:

  Logical. If TRUE, the function will return the combined data frame.
  Default is FALSE, meaning it will only write the file and not return
  the data.

## Value

If `ret` is TRUE, a data frame containing the combined data from all
input files. If `ret` is FALSE, the function writes the combined data to
a TSV file named `ipr_combined.tsv` in the specified directory and
returns NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
combineIPR <- combine_ipr("path/to/ipr/files", ret = TRUE)
} # }
```
