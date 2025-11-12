# createRepresentativeAccNum

Function to generate a vector of one Accession number per distinct
observation from 'reduced' column

Function to generate a vector of one Accession number per distinct
observation from 'reduced' column

## Usage

``` r
createRepresentativeAccNum(
  prot_data,
  reduced = "Lineage",
  accnum_col = "AccNum"
)

createRepresentativeAccNum(
  prot_data,
  reduced = "Lineage",
  accnum_col = "AccNum"
)
```

## Arguments

- prot_data:

  Data frame containing Accession Numbers

- reduced:

  Column from prot_data from which distinct observations will be
  generated from. One accession number will be assigned for each of
  these observations

- accnum_col:

  Column from prot_data that contains Accession Numbers

## Value

A character vector containing one Accession number per distinct
observation from the specified reduced column.

A character vector containing representative accession numbers, one for
each distinct observation in the specified 'reduced' column.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
createRepresentativeAccNum(prot)
} # }
if (FALSE) { # \dontrun{
# Example usage with a data frame called `protein_data`
createRepresentativeAccNum <- RepresentativeAccNums(prot_data = protein_data,
                                                    reduced = "Lineage",
                                                    accnum_col = "AccNum")
representative_accessions
} # }
```
