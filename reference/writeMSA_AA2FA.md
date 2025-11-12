# writeMSA_AA2FA

MsaAAMultipleAlignment Objects are generated from calls to
msaClustalOmega and msaMuscle from the 'msa' package

Write MsaAAMultpleAlignment Objects as aligned fasta sequence
MsaAAMultipleAlignment Objects are generated from calls to
msaClustalOmega and msaMuscle from the 'msa' package

## Usage

``` r
writeMSA_AA2FA(alignment, outpath)

writeMSA_AA2FA(alignment, outpath)
```

## Arguments

- alignment:

  MsaAAMultipleAlignment object to be written as a fasta

- outpath:

  Where the resulting FASTA file should be written to

## Value

Character string representing the content of the written FASTA file.

Character string of the FASTA content that was written to the file.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
writeMSA_AA2FA("my_sequences.fasta", outpath = "aligned_output.fasta")
} # }
if (FALSE) { # \dontrun{
# Example usage
alignment <- alignFasta("path/to/sequences.fasta")
writeMSA_AA2FA(alignment, "path/to/aligned_sequences.fasta")
} # }
```
