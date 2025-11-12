# runIPRScan

Run InterProScan on a given FASTA file and save the results to an output
file.

## Usage

``` r
runIPRScan(filepath_fasta, filepath_out, appl = c("Pfam", "Gene3D"))
```

## Arguments

- filepath_fasta:

  A string representing the path to the input FASTA file.

- filepath_out:

  A string representing the base path for the output file.

- appl:

  A character vector specifying the InterProScan applications to use
  (e.g., "Pfam", "Gene3D"). Default is `c("Pfam", "Gene3D")`.

## Value

A data frame containing the results from the InterProScan output TSV
file.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- runIPRScan(
    filepath_fasta = "path/to/your_fasta_file.fasta",
    filepath_out = "path/to/output_file",
    appl = c("Pfam", "Gene3D")
)
results
} # }
```
