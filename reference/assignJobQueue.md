# assignJobQueue

Decision function to assign job queue

## Usage

``` r
assignJobQueue(t_sec_estimate, t_cutoff = 21600)
```

## Arguments

- t_sec_estimate:

  estimated number of seconds a job will process (from
  calculateEstimatedWallTimeFromOpts ())

- t_long:

  threshold value that defines the lower bound for assigning a job to
  the "long queue"

## Value

a string of "short" or "long"

example: calculateEstimatedWallTimeFromOpts (c("homology_search",
"domain_architecture"), 3) \|\> assignJobQueue()
