# Download the combined assembly summaries of genbank and refseq

Download the combined assembly summaries of genbank and refseq

## Usage

``` r
combineFiles(
  inpath = c("../molevol_data/project_data/phage_defense/"),
  pattern = "*full_analysis.tsv",
  delim = "\t",
  skip = 0,
  col_names = T
)
```

## Arguments

- inpath:

  Character. The master directory path where the files reside. The
  search is recursive (i.e., it will look in subdirectories as well).

- pattern:

  Character. A search pattern to identify files to be combined. Default
  is "\*full_analysis.tsv".

- delim:

  Character. The delimiter used in the input files. Default is tab ("").

- skip:

  Integer. The number of lines to skip at the beginning of each file.
  Default is 0.

- col_names:

  Logical or character vector. If TRUE, the first row of each file is
  treated as column names. Alternatively, a character vector can be
  provided to specify custom column names.

## Value

A data frame containing the combined contents of all matched files. Each
row will include a new column "ByFile" indicating the source file of the
data.

## Author

Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
combined_data <- combineFiles(inpath = "../molevol_data/project_data/phage_defense/")
} # }
```
