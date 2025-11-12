# acc2FA

converts protein accession numbers to a fasta format. Resulting fasta
file is written to the outpath.

converts protein accession numbers to a fasta format. Resulting fasta
file is written to the outpath.

acc2FA converts protein accession numbers to a fasta format. Resulting
fasta file is written to the outpath.

## Usage

``` r
acc2FA(accessions, outpath, plan = "sequential")

acc2FA(accessions, outpath, plan = "sequential")
```

## Arguments

- accessions:

  Character vector containing protein accession numbers to generate
  fasta sequences for. Function may not work for vectors of length \>
  10,000

- outpath:

  [str](https://rdrr.io/r/utils/str.html). Location where fasta file
  should be written to.

- plan:

  Character. The plan to use for processing. Default is "sequential".

## Value

A logical value indicating whether the retrieval and conversion were
successful. Returns `TRUE` if successful and `FALSE` otherwise.

A Fasta file is written to the specified `outpath`.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
acc2FA(accessions = c("ACU53894.1", "APJ14606.1", "ABK37082.1"),
outpath = "my_proteins.fasta")
Entrez:accessions <- rep("ANY95992.1", 201) |> acc2FA(outpath = "entrez.fa")
EBI:accessions <- c("P12345", "Q9UHC1",
"O15530", "Q14624", "P0DTD1") |> acc2FA(outpath = "ebi.fa")
} # }
if (FALSE) { # \dontrun{
acc2FA(accessions = c("ACU53894.1", "APJ14606.1", "ABK37082.1"),
outpath = "my_proteins.fasta")
Entrez:accessions <- rep("ANY95992.1", 201) |> acc2FA(outpath = "entrez.fa")
EBI:accessions <- c("P12345", "Q9UHC1", "O15530", "Q14624", "P0DTD1") |>
acc2FA(outpath = "ebi.fa")
} # }
```
