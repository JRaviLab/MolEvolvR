# writeProcessRuntime2YML

Compute median process runtimes, then write a YAML list of the processes
and their median runtimes in seconds to the path specified by
'filepath'.

The default value of filepath is the value of the env var
MOLEVOLVR_PROC_WEIGHTS, which getProcessRuntimeWeights() also uses as
its default read location.

## Usage

``` r
writeProcessRuntime2YML(dir_job_results, filepath = NULL)
```

## Arguments

- dir_job_results:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  path to MolEvolvR job_results directory

- filepath:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  path to save YAML file; if NULL, uses
  ./molevol_scripts/log_data/job_proc_weights.yml

## Examples

``` r
if (FALSE) { # \dontrun{
writeProcessRuntime2YML(
    "/data/scratch/janani/molevolvr_out/",
    "/data/scratch/janani/molevolvr_out/log_tbl.yml"
)
} # }
```
