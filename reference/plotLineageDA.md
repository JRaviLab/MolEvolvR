# Lineage Plot: Heatmap of Domains/DAs/GCs vs Lineages

Lineage plot for Domains, Domain Architectures and Genomic Contexts.
Heatmap.

## Usage

``` r
plotLineageDA(
  query_data = "prot",
  colname = "DomArch",
  cutoff = 90,
  RowsCutoff = FALSE,
  color = "default"
)
```

## Arguments

- query_data:

  Data frame of protein homologs with the usual 11 columns + additional
  word columns (0/1 format). Default is prot (variable w/ protein data).

- colname:

  Column name from query_data: "DomArch.norep", "GenContext.norep",
  "DomArch.PFAM.norep" or "DomArch.LADB.norep". Default is
  "DomArch.norep".

- cutoff:

  Numeric. Cutoff for word frequency. Default is 90.

- RowsCutoff:

  Boolean. If TRUE, applies a row cutoff to remove data rows based on a
  certain condition. Default is FALSE.

- color:

  Color for the heatmap. One of six options: "default", "magma",
  "inferno", "plasma", "viridis", or "cividis"

## Value

A LineageDA plot object.

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
plotLineageDA(toast_rack_data, 10, "DomArch.norep", "da2doms")
} # }
```
