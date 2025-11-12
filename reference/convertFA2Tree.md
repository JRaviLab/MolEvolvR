# convertFA2Tree

convertFA2Tree

## Usage

``` r
convertFA2Tree(
  fa_path = here("data/alns/pspa_snf7.fa"),
  tre_path = here("data/alns/pspa_snf7.tre"),
  fasttree_path = here("src/FastTree")
)
```

## Arguments

- fa_path:

  Path to the input FASTA alignment file (.fa). Default is the path to
  "data/alns/pspa_snf7.fa".

- tre_path:

  Path to the output file where the generated tree (.tre) will be saved.
  Default is the path to "data/alns/pspa_snf7.tre".

- fasttree_path:

  Path to the FastTree executable, which is used to generate the
  phylogenetic tree. Default is "src/FastTree".

## Value

No return value. The function generates a tree file (.tre) from the
input FASTA file.

## Examples
