# calculateProcessRuntime

Scrape MolEvolvR logs and calculate median processes

## Usage

``` r
calculateProcessRuntime(dir_job_results)
```

## Arguments

- dir_job_results:

  [rlang::chr](https://rlang.r-lib.org/reference/vector-construction.html)
  path to MolEvolvR job_results directory

## Value

[list](https://rdrr.io/r/base/list.html) names: processes; values:
median runtime (seconds)

see molevol_scripts/R/metrics.R for info on functions called here

examples:

1.  

dir_job_results \<- "/data/scratch/janani/molevolvr_out"
list_proc_medians \<- calculateProcessRuntime(dir_job_results)

1.  from outside container environment common_root \<-
    "/data/molevolvr_transfer/molevolvr_dev" dir_job_results \<-
    "/data/molevolvr_transfer/molevolvr_dev/job_results"
    list_proc_medians \<- calculateProcessRuntime(dir_job_results)
