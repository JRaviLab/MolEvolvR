# calculateEstimatedWallTimeFromOpts

Given MolEvolvR advanced options and number of inputs, calculate the
total estimated walltime for the job

## Usage

``` r
calculateEstimatedWallTimeFromOpts(
  advanced_opts,
  n_inputs = 1L,
  n_hits = NULL,
  verbose = FALSE
)
```

## Arguments

- advanced_opts:

  character vector of MolEvolvR advanced options (see mapOption2Process
  for the options)

- n_inputs:

  total number of input proteins

## Value

total estimated number of seconds a job will process (walltime)

example: calculateEstimatedWallTimeFromOpts (c("homology_search",
"domain_architecture"), n_inputs = 3, n_hits = 50L)
