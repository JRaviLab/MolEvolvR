# Multiple Sequence Alignment

Generates a multiple sequence alignment from a fasta file

createMSA_PDF is a function that reads a fasta file and generates a
multiple sequence alignment as a pdf

## Usage

``` r
createMSA_PDF(
  fasta_path,
  out_path = NULL,
  lowerbound = NULL,
  upperbound = NULL
)
```

## Arguments

- fasta_path:

  Character. The path location of the fasta file to be read.

- out_path:

  Character. The path location of the output pdf to write. Default is
  NULL. If value is NULL, the pdf is written to the same directory as
  the fasta file.

- lowerbound:

  Numeric. The column that determines the starting location of the MSA.
  Default is NULL. If value is NULL, the entire multiple sequence
  alignment is printed.

- upperbound:

  Numeric. The column that determines the ending location of the MSA.
  Default is NULL. If value is NULL, the entire multiple sequence
  alignment is printed.

## Value

A PDF file containing the multiple sequence alignment.

## Examples

``` r
if (FALSE) { # \dontrun{
createMSA_PDF(fasta_path = "path/to/your/file.fasta", 
        out_path = "path/to/output/alignment.pdf", 
        lowerbound = 10, 
        upperbound = 200)
} # }
```
