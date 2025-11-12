# proteinAcc2TaxID_old

Perform elink to go from protein database to taxonomy database and write
the resulting file of taxid and lineage to out_path

## Usage

``` r
proteinAcc2TaxID_old(accessions, out_path, plan = "multicore")
```

## Arguments

- accessions:

  A character vector containing the accession numbers to query in the
  protein database.

- out_path:

  A string specifying the path where the results of the query will be
  written. If set to NULL, a temporary directory will be used.

- plan:

  A character string that specifies the execution plan for parallel
  processing. The default is "multicore".

## Value

This function does not return a value. It writes the results to the
specified output path.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
accessions <- c("ABC123", "XYZ456", "LMN789")
proteinAcc2TaxID_old(accessions, out_path = "/path/to/output")
} # }
```
