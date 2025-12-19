# convertAlignment2FA

convertAlignment2FA

## Usage

``` r
convertAlignment2FA_legacy(
  aln_file = "",
  lin_file = "data/rawdata_tsv/all_semiclean.txt",
  fa_outpath = "",
  reduced = FALSE
)
```

## Arguments

- aln_file:

  Character. Path to file. Input tab-delimited file + alignment file
  accnum & alignment. Default is 'pspa_snf7.aln'

- lin_file:

  Character. Path to file. Protein file with accession + number to
  lineage mapping. Default is 'pspa.txt'

- fa_outpath:

  Character. Path to the written fasta file. Default is 'NULL'

- reduced:

  Boolean. If TRUE, the fasta file will contain only one sequence per
  lineage. Default is 'FALSE'

## Value

Character string containing the Fasta formatted sequences. If
`fa_outpath` is specified, the function also writes the sequences to the
Fasta file.

## Details

The alignment file would need two columns: 1. accession + number and 2.
alignment. The protein homolog accession to lineage mapping + file
should have

## Note

Please refer to the source code if you have alternate + file formats
and/or column names.

## Author

Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
convertAlignment2FA("pspa_snf7.aln", "pspa.txt")
} # }
```
