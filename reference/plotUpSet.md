# UpSet Plot

UpSet plot for Domain Architectures vs Domains and Genomic Contexts vs
Domain Architectures.

## Usage

``` r
plotUpSet(
  query_data = "toast_rack.sub",
  colname = "DomArch",
  cutoff = 90,
  RowsCutoff = FALSE,
  text.scale = 1.5,
  point.size = 2.2,
  line.size = 0.8
)
```

## Arguments

- query_data:

  Data frame of protein homologs with the usual 11 columns + additional
  word columns (0/1 format). Default is toast_rack.sub

- colname:

  Column name from query_data: "DomArch.norep", "GenContext.norep",
  "DomArch.PFAM.norep" or "DomArch.LADB.norep". Default is
  "DomArch.norep".

- cutoff:

  Numeric. Cutoff for word frequency. Default is 90.

- RowsCutoff:

  Boolean. If TRUE, applies a row cutoff to remove data rows based on a
  certain condition. Default is FALSE.

- text.scale:

  Allows scaling of axis title, tick lables, and numbers above the
  intersection size bars. text.scale can either take a universal scale
  in the form of an integer, or a vector of specific scales in the
  format: c(intersection size title, intersection size tick labels, set
  size title, set size tick labels, set names, numbers above bars)

- point.size:

  Numeric. Sets the size of points in the UpSet plot. Default is 2.2.

- line.size:

  Numeric. Sets the line width in the UpSet plot. Default is 0.8.

## Value

An UpSet plot object. The plot visualizes intersections of sets based on
the provided colname in query_data.

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
plotUpSet(pspa.sub, 10, "da2doms")
} # }
```
