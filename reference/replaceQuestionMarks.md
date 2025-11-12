# Replace QMs

Replace consecutive '?' separated by '-\>', '\<-' or '\|\|' with 'X(s)'
Replace '?' with 'X'

## Usage

``` r
replaceQuestionMarks(prot, by_column = "GenContext")
```

## Arguments

- prot:

  DataTable to operate on

- by_column:

  Column to operate on

## Value

The original data frame with the specified column updated. All
consecutive '?' characters will be replaced with 'X(s)', and individual
'?' characters will be replaced with 'X'.

## Examples

``` r
if (FALSE) { # \dontrun{
replaceQuestionMarks()
} # }
```
