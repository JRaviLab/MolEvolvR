# mapAdvOption2Process

Use MolEvolvR advanced options to get associated processes

## Usage

``` r
mapAdvOption2Process(advanced_opts)
```

## Arguments

- advanced_opts:

  character vector of MolEvolvR advanced options

## Value

character vector of process names that will execute given the advanced
options

example: advanced_opts \<- c("homology_search", "domain_architecture")
procs \<- mapAdvOption2Process(advanced_opts)
