# efetchIPG

Perform efetch on the ipg database and write the results to out_path

Perform efetch on the ipg database and write the results to out_path

## Usage

``` r
efetchIPG(accessions, out_path, plan = "multicore")

efetchIPG(accessions, out_path, plan = "multicore")
```

## Arguments

- accessions:

  Character vector containing the accession numbers to query on the ipg
  database

- out_path:

  Path to write the efetch results to

- plan:

  Character. Specifies the execution plan for parallel processing.
  Default is "multicore".

- accnums:

  Character vector containing the accession numbers to query on the ipg
  database

## Value

No return value. The function writes the fetched results to `out_path`.

The function does not return a value but writes the efetch results
directly to the specified `out_path`.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
efetchIPG()
} # }
if (FALSE) { # \dontrun{
efetchIPG(
  accessions = c("P12345", "Q67890", "A12345"),
  out_path = "path/to/efetch_results.xml"
)
} # }
```
