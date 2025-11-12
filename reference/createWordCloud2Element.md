# Wordclouds for the predominant domains, domain architectures.

Wordclouds for the predominant domains (from DAs) and DAs (from GC)

## Usage

``` r
createWordCloud2Element(
  query_data = "prot",
  colname = "DomArch",
  cutoff = 70,
  UsingRowsCutoff = FALSE
)
```

## Arguments

- query_data:

  Data frame of protein homologs with the usual 11 columns + additional
  word columns (0/1 format). Default is "prot".

- colname:

  Character. The name of the column in `query_data` to generate the word
  cloud from. Default is "DomArch".

- cutoff:

  Numeric. The cutoff value for filtering elements based on their
  frequency. Default is 70.

- UsingRowsCutoff:

  Logical. Whether to use a row-based cutoff instead of a frequency
  cutoff. Default is FALSE.

## Value

A word cloud plot showing the frequency of elements from the selected
column.

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
createWordCloudElement(prot, "da2doms", 10)
} # }
```
