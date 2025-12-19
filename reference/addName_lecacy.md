# addName

This function adds a new 'Name' column that is comprised of components
from Kingdom, Phylum, Genus, and species, as well as the accession

## Usage

``` r
addName_lecacy(
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

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
addName(data_frame)
} # }
```
