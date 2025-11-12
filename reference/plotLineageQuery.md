# Lineage Plot: Heatmap of Queries vs Lineages

Lineage plot for queries. Heatmap.

## Usage

``` r
plotLineageQuery(
  query_data = all,
  queries,
  colname = "ClustName",
  cutoff,
  color = "default"
)
```

## Arguments

- query_data:

  Data frame of protein homologs with the usual 11 columns + additional
  word columns (0/1 format). Default is prot (variable w/ protein data).

- queries:

  Character Vector containing the queries that will be used for the
  categories.

- colname:

  Character. The column used for filtering based on the `queries`.
  Default is "ClustName".

- cutoff:

  Numeric. The cutoff value for filtering rows based on their total
  count. Rows with values below this cutoff are excluded.

- color:

  Character. Defines the color palette used for the heatmap. Default is
  a red gradient.

## Value

A ggplot object representing a heatmap (tile plot) showing the
relationship between queries and lineages, with the intensity of color
representing the count of matching records.

## Note

Please refer to the source code if you have alternate file formats
and/or column names.

## Author

Janani Ravi, Samuel Chen

## Examples

``` r
if (FALSE) { # \dontrun{
plotLineageQuery(prot, c("PspA", "PspB", "PspC", "PspM", "PspN"), 95)
} # }
```
