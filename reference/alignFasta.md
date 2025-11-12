# alignFasta

Perform a Multiple Sequence Alignment on a FASTA file.

Perform a Multiple Sequence Alignment on a FASTA file.

## Usage

``` r
alignFasta(fasta_file, tool = "Muscle", outpath = NULL)

alignFasta(fasta_file, tool = "Muscle", outpath = NULL)
```

## Arguments

- fasta_file:

  Path to the FASTA file to be aligned

- tool:

  Type of alignment tool to use. One of three options: "Muscle",
  "ClustalO", or "ClustalW"

- outpath:

  Path to write the resulting alignment to as a FASTA file. If NULL, no
  file is written

## Value

aligned fasta sequence as a MsaAAMultipleAlignment object

aligned fasta sequence as a MsaAAMultipleAlignment object

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
aligned_sequences <- alignFasta("my_sequences.fasta",
tool = "Muscle", outpath = "aligned_output.fasta")
} # }
if (FALSE) { # \dontrun{
# Example usage
aligned_sequences <- alignFasta("path/to/sequences.fasta",
tool = "ClustalO", outpath = "path/to/aligned_sequences.fasta")
aligned_sequences
} # }
```
