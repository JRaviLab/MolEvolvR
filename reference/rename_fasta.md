# Rename the labels of fasta files

Rename the labels of fasta files

## Usage

``` r
rename_fasta(fa_path, outpath, replacement_function = .data$map_acc2name, ...)
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

Character vector containing the modified lines of the Fasta file.

## Examples

``` r
if (FALSE) { # \dontrun{
rename_fasta("input.fasta", "output.fasta",
replacement_function = map_acc2name, acc2name = acc2name_table)
} # }
```
