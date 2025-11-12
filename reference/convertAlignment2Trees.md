# convertAlignment2Trees

Generate Trees for ALL fasta files in "data/alns"

## Usage

``` r
convertAlignment2Trees(aln_path = here("data/alns/"))
```

## Arguments

- aln_path:

  Path to the directory containing all the alignment FASTA files (.fa)
  for which trees will be generated. Default is "data/alns/".

## Value

No return value. The function generates tree files (.tre) for each
alignment file in the specified directory.

## Examples

``` r
if (FALSE) { # \dontrun{
generate_trees(here("data/alns/"))
} # }
```
