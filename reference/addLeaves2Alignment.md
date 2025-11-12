# addLeaves2Alignment

Adding Leaves to an alignment file w/ accessions Genomic Contexts vs
Domain Architectures.

Adding Leaves to an alignment file w/ accessions Genomic Contexts vs
Domain Architectures.

## Usage

``` r
addLeaves2Alignment(
  aln_file = "",
  lin_file = "data/rawdata_tsv/all_semiclean.txt",
  reduced = FALSE
)

addLeaves2Alignment(
  aln_file = "",
  lin_file = "data/rawdata_tsv/all_semiclean.txt",
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

- reduced:

  Boolean. If TRUE, a reduced data frame will be generated with only one
  sequence per lineage. Default is FALSE.

## Value

A data frame containing the enriched alignment data with lineage
information.

A data frame containing the combined alignment and lineage information.

## Details

The alignment file would need two columns: 1. accession + number and 2.
alignment. The protein homolog accession to lineage mapping + file
should have

The alignment file would need two columns: 1. accession + number and 2.
alignment. The protein homolog accession to lineage mapping + file
should have

## Note

Please refer to the source code if you have alternate + file formats
and/or column names.

Please refer to the source code if you have alternate + file formats
and/or column names.

## Author

Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
addLeaves2Alignment("pspa_snf7.aln", "pspa.txt")
} # }
if (FALSE) { # \dontrun{
addLeaves2Alignment("pspa_snf7.aln", "pspa.txt")
} # }
```
