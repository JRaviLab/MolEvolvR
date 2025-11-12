# Function to generate MSA using kalign

Function to generate MSA using kalign

## Usage

``` r
createMSA_Kalign(fa_file = "", outfile = "")
```

## Arguments

- fa_file:

  Character. The path to the input FASTA file containing protein
  sequences.

- outfile:

  Character. The path to the output file where the alignment will be
  saved.

## Value

A list containing the alignment object and the output file path.

## Examples

``` r
if (FALSE) { # \dontrun{
createMSA_Kalign(fa_file = "path/to/sequences.fasta", 
                 outfile = "path/to/alignment.txt")
} # }
```
