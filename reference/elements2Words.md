# elements2Words

Break string ELEMENTS into WORDS for domain architecture (DA) and
genomic context (GC)

## Usage

``` r
elements2Words(prot, column = "DomArch", conversion_type = "da2doms")
```

## Arguments

- prot:

  A dataframe containing the dataset to analyze. The specified `column`
  contains the string elements to be processed.

- column:

  A character string specifying the name of the column to analyze.
  Default is "DomArch".

- conversion_type:

  A character string specifying the type of conversion. Two options are
  available:

  `da2doms`

  :   Convert domain architectures into individual domains by replacing
      `+` symbols with spaces.

  `gc2da`

  :   Convert genomic context into domain architectures by replacing
      directional symbols (`<-`, `->`, and `|`) with spaces.

## Value

A single string where elements are delimited by spaces. The function
performs necessary substitutions based on the `conversion_type` and
cleans up extraneous characters like newlines, tabs, and multiple
spaces.

## Examples

``` r
if (FALSE) { # \dontrun{
tibble::tibble(DomArch = c("aaa+bbb",
"a+b", "b+c", "b-c")) |> elements2Words()
} # }
```
