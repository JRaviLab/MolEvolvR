# Lineage Plot for top neighbors

Lineage plot for top neighbors obtained from DAs of Genomic Contexts.

## Usage

``` r
plotLineageNeighbors(
  query_data = "prot",
  query = "pspa",
  colname = "GenContext.norep"
)
```

## Arguments

- query_data:

  Data frame of protein homologs with the usual 11 columns + additional
  word columns (0/1 format). Default is pspa_data.

- query:

  Name of query protein/domain. Default is "pspa".

- colname:

  Column name from query_data. Default is "GenContext.norep".

## Value

A ggplot object representing a heatmap (tile plot) of lineage versus the
top neighboring domain architectures, with color intensity representing
the frequency of occurrences.

## Details

For "da2doms" you would need the file DA.doms.wc as well as the column
query_data\$DomArch.norep

For "gc2da", you would need the file GC.DA.wc as well as the column
query_data\$GenContext.norep

## Note

Please refer to the source code if you have alternate file formats
and/or column names.

## Author

Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
plotLineageNeighbors(pspa_data, pspa, "GenContext.norep", "da2doms")
} # }
```
