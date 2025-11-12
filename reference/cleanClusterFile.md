# Clean Cluster File

Reads and cleans a cluster file

This function reads a space-separated cluster file and converts it to a
cleaned up data frame.

## Usage

``` r
cleanClusterFile(path, writepath = NULL, query)
```

## Arguments

- path:

  A character to the path of the cluster file to be cleaned

- writepath:

  A character designating where the tsv file of the cleaned cluster file
  will be written to. If value is NULL no file is written. Default NULL

- query:

  A character identifying the query of the file.

## Value

The cleaned up cluster data frame is returned and a tsv file is written
if the "writepath" parameter is used.

## Examples

``` r
if (FALSE) { # \dontrun{
cleanClusterFile("data/pspa.op_ins_cls", writepath = NULL, query = "pspa")
} # }
```
