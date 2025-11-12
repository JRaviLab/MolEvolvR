# Combining full_analysis files

Combining full_analysis files

## Usage

``` r
combineFullAnalysis(inpath, ret = FALSE)
```

## Arguments

- inpath:

  Character. The path to the directory containing the
  `.full_analysis.tsv` files to be combined.

- ret:

  Logical. If TRUE, the function will return the combined data frame.
  Default is FALSE, meaning it will only write the file and not return
  the data.

## Value

If `ret` is TRUE, a data frame containing the combined data from all
input files. If `ret` is FALSE, the function writes the combined data to
a TSV file named `cln_combined.tsv` in the specified directory and
returns NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
combined_data <- combineFullAnalysis("path/to/full_analysis/files", ret = TRUE)
} # }
```
