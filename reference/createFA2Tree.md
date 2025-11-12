# createFA2Tree

Generating phylogenetic tree from alignment file '.fa'

## Usage

``` r
createFA2Tree(
  fa_file = "data/alns/pspa_snf7.fa",
  out_file = "data/alns/pspa_snf7.tre"
)
```

## Arguments

- fa_file:

  Character. Path to the alignment FASTA file (.fa) from which the
  phylogenetic tree will be generated. Default is 'pspa_snf7.fa'.

- out_file:

  Path to the output file where the generated tree (.tre) will be saved.
  Default is "data/alns/pspa_snf7.tre".

## Value

No return value. The function generates a phylogenetic tree file (.tre)
based on different approaches like Neighbor Joining, UPGMA, and Maximum
Likelihood.

## Details

The alignment file would need two columns: 1. accession + number and 2.
alignment. The protein homolog accession to lineage mapping + file
should have

## Note

Please refer to the source code if you have alternate + file formats
and/or column names.

## Author

Janani Ravi, MolEcologist

## Examples

``` r
if (FALSE) { # \dontrun{
generate_aln2tree("pspa_snf7.fa")
} # }
```
