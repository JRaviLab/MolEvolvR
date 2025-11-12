# addTaxID

addTaxID

## Usage

``` r
addTaxID(data, acc_col = "AccNum", version = T)
```

## Arguments

- data:

  A data frame or data table containing protein accession numbers.

- acc_col:

  A string specifying the column name in `data` that contains the
  accession numbers. Defaults to "AccNum".

- version:

  A logical indicating whether to remove the last two characters from
  the accession numbers for TaxID retrieval. Defaults to TRUE.

## Value

A data table that includes the original data along with a new column
containing the corresponding TaxIDs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a sample data table with accession numbers
sample_data <- data.table(AccNum = c("ABC123.1", "XYZ456.1", "LMN789.2"))
enriched_data <- addTaxID(sample_data, acc_col = "AccNum", version = TRUE)
enriched_data
} # }
```
