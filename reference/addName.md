# addName

This function adds a new 'Name' column that is comprised of components
from Kingdom, Phylum, Genus, and species, as well as the accession

This function adds a new 'Name' column that is comprised of components
from Kingdom, Phylum, Genus, and species, as well as the accession

## Usage

``` r
addName(
  data,
  accnum_col = "AccNum",
  spec_col = "Species",
  lin_col = "Lineage",
  lin_sep = ">",
  out_col = "Name"
)

addName(
  data,
  accnum_col = "AccNum",
  spec_col = "Species",
  lin_col = "Lineage",
  lin_sep = ">",
  out_col = "Name"
)
```

## Arguments

- data:

  Data to add name column to

- accnum_col:

  Column containing accession numbers

- spec_col:

  Column containing species

- lin_col:

  Column containing lineage

- lin_sep:

  Character separating lineage levels

- out_col:

  Column that contains the new 'Name' derived from Species, Lineage, and
  AccNum info

## Value

Original data with a 'Name' column

Original data with a 'Name' column

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
# Example usage of the addName function
data <- data.frame(
  AccNum = c("ACC123", "ACC456"),
  Species = c("Homo sapiens", "Mus musculus"),
  Lineage = c("Eukaryota>Chordata", "Eukaryota>Chordata")
)
enriched_data <- addName(data)
enriched_data
#>   AccNum      Species            Lineage                     Name
#> 1 ACC123 Homo sapiens Eukaryota>Chordata  EChorda_Hsapiens_ACC123
#> 2 ACC456 Mus musculus Eukaryota>Chordata EChorda_Mmusculus_ACC456
if (FALSE) { # \dontrun{
addName(data_frame)
} # }
```
