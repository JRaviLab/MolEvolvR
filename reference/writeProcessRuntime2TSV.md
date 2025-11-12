# writeProcessRuntime2TSV

Write a table of 2 columns: 1) process and 2) median seconds

## Usage

``` r
writeProcessRuntime2TSV(dir_job_results, filepath)
```

## Arguments

- dir_job_results:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  path to MolEvolvR job_results

- filepath:

  path to save tsv file

## Value

[tbl_df](https://dplyr.tidyverse.org/reference/tbl_df.html) 2
columns: 1) process and 2) median seconds

example: writeProcessRuntime2TSV( "/data/scratch/janani/molevolvr_out/",
"/data/scratch/janani/molevolvr_out/log_tbl.tsv" )
