# getTopAccByLinDomArch

Group by lineage + DA then take top 20

## Usage

``` r
getTopAccByLinDomArch(
  infile_full,
  DA_col = "DomArch.Pfam",
  lin_col = "Lineage_short",
  n = 20,
  query
)
```

## Arguments

- infile_full:

  A data frame containing the full dataset with lineage and domain
  architecture information.

- DA_col:

  A string representing the name of the domain architecture column.
  Default is "DomArch.Pfam".

- lin_col:

  A string representing the name of the lineage column. Default is
  "Lineage_short".

- n:

  An integer specifying the number of top accession numbers to return.
  Default is 20.

- query:

  A string for filtering a specific query name. If it is not "All", only
  the data matching this query will be processed.

## Value

A vector of the top N accession numbers (`AccNum`) based on counts
grouped by lineage and domain architecture.

## Examples

``` r
if (FALSE) { # \dontrun{
top_accessions <- getTopAccByLinDomArch(infile_full = my_data,
DA_col = "DomArch.Pfam", lin_col = "Lineage_short",
n = 20, query = "specific_query_name")
} # }
```
