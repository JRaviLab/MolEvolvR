# Run DELTABLAST to find homologs for proteins of interest

This function executes a Delta-BLAST search using the specified
parameters and database. It sets the BLAST database path, runs the
Delta-BLAST command with the given query, and outputs the results.

## Usage

``` r
runDeltaBlast(
  runDeltaBlast,
  db_search_path,
  db = "refseq",
  query,
  evalue = "1e-5",
  out,
  num_alignments,
  num_threads = 1
)
```

## Arguments

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

- num_alignments:

  Number of alignments to report.

- num_threads:

  Number of threads to use for the search (default is 1).

- deltablast_path:

  Path to the Delta-BLAST executable.

## Value

This function does not return a value; it outputs results to the
specified file.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
runDeltaBlast(runDeltaBlast, db_search_path)
} # }
```
