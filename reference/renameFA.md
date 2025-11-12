# Rename the labels of fasta files

Rename the labels of fasta files

## Usage

``` r
renameFA(fa_path, outpath, replacement_function = mapAcc2Name, ...)
```

## Arguments

- fa_path:

  Path to fasta file

- outpath:

  Path to write altered fasta file to

- replacement_function:

  Function to apply to lines starting with '\>'

- ...:

  Additional arguments to pass to replacement_function

## Value

A character vector of the modified lines in the FASTA file.

## Examples

``` r
if (FALSE) { # \dontrun{
renameFA("path/to/input.fasta",
"path/to/output.fasta", mapAcc2Name, acc2name)
} # }
```
