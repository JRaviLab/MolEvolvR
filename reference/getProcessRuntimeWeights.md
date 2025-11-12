# getProcessRuntimeWeights

Quickly get the runtime weights for MolEvolvR backend processes

## Usage

``` r
getProcessRuntimeWeights(medians_yml_path = NULL)
```

## Arguments

- dir_job_results:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  path to MolEvolvR job_results directory

## Value

[list](https://rdrr.io/r/base/list.html) names: processes; values:
median runtime (seconds)

example: writeProcessRuntime2YML()
