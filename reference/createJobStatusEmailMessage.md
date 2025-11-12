# Produces a mail message that can be sent to a user when their job is accepted. Used by the sendJobStatusEmail() method.

Produces a mail message that can be sent to a user when their job is
accepted. Used by the sendJobStatusEmail() method.

## Usage

``` r
createJobStatusEmailMessage(
  job_dir,
  pin_id,
  job_results_url,
  event_type,
  context
)
```

## Arguments

- job_dir:

  the directory where the job's arguments are stored, in job_args.yml

- pin_id:

  the unique identifier for the job

- job_results_url:

  the URL where the user can check the status of their job

- event_type:

  either 'start' or 'end', returns the corresponding email for the given
  type

- context:

  a list of additional values, e.g. job runtime info, that can be used
  in the template emails

## Value

the result of the sendmailR::sendmail() call
