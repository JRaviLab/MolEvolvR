# Run RPSBLAST to generate domain architectures for proteins of interest

This function executes an RPS-BLAST search to generate domain
architectures for specified proteins. It sets the BLAST database path,
runs the RPS-BLAST command with the provided query, and outputs the
results.

## Usage

``` r
runRPSBlast(
  rpsblast_path,
  db_search_path,
  db = "refseq",
  query,
  evalue = "1e-5",
  out,
  num_threads = 1
)
```

## Arguments

- rpsblast_path:

  Path to the RPS-BLAST executable.

- db_search_path:

  Path to the BLAST databases.

- db:

  Name of the BLAST database to search against (default is "refseq").

- query:

  Path to the input query file.

- evalue:

  E-value threshold for reporting matches (default is "1e-5").

- out:

  Path to the output file where results will be saved.

- num_threads:

  Number of threads to use for the search (default is 1).

## Value

This function does not return a value; it outputs results to the
specified file.

## Examples

``` r
if (FALSE) { # \dontrun{
runRSPBlast(rpsblast_path, db_search_path, query, out)
} # }
```
