# generateAllAlignments2FA

Adding Leaves to an alignment file w/ accessions

Adding Leaves to all alignment files w/ accessions & DAs?

## Usage

``` r
generateAllAlignments2FA_legacy(
  aln_path = here("data/rawdata_aln/"),
  fa_outpath = here("data/alns/"),
  lin_file = here("data/rawdata_tsv/all_semiclean.txt"),
  reduced = F
)
```

## Arguments

- aln_path:

  Character. Path to alignment files. Default is
  'here("data/rawdata_aln/")'

- fa_outpath:

  Character. Path to file. Master protein file with AccNum & lineages.
  Default is 'here("data/rawdata_tsv/all_semiclean.txt")'

- lin_file:

  Character. Path to the written fasta file. Default is
  'here("data/alns/")'.

- reduced:

  Boolean. If TRUE, the fasta file will contain only one sequence per
  lineage. Default is 'FALSE'.

## Value

NULL. The function saves the output FASTA files to the specified
directory.

## Details

The alignment files would need two columns separated by spaces:

1.  AccNum and 2. alignment. The protein homolog file should have
    AccNum, Species, Lineages.

## Note

Please refer to the source code if you have alternate + file formats
and/or column names.

## Examples

``` r
if (FALSE) { # \dontrun{
generateAllAlignments2FA()
} # }
```
