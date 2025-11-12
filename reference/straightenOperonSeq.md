# straightenOperonSeq: Reverse Equalities in Genomic Context

This function processes the genomic context strings (GenContext) and
reverses directional signs based on the presence of an equal sign ("=").

## Usage

``` r
straightenOperonSeq(prot)
```

## Arguments

- prot:

  [vector](https://rdrr.io/r/base/vector.html) A vector of genomic
  context strings to be processed.

## Value

[vector](https://rdrr.io/r/base/vector.html) A vector of the same length
as the input, where each genomic element is annotated with either a
forward ("-\>") or reverse ("\<-") direction, depending on its position
relative to the "=" symbols.

## Examples

``` r
# Example input: Genomic context with directional symbols and an asterisk
genomic_context <- c("A", "B", "*", "C", "D", "=", "E", "F")
straightenOperonSeq(genomic_context)
#> [1] "A->" "B->" "*->" "C->" "D->" "="   "<-E" "<-F"

# Output: "A->", "B->", "*", "<-C", "<-D", "=", "E->", "F->"
```
