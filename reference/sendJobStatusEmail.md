# Sends a "job accepted" email to a user when their job is accepted, including details about the job submission and how to check its status.

Sends a "job accepted" email to a user when their job is accepted,
including details about the job submission and how to check its status.

## Usage

``` r
sendJobStatusEmail(notify_email, job_dir, pin_id, event_type, context = NULL)
```

## Arguments

- notify_email:

  the email address to send the notification to

- job_dir:

  the directory where the job's arguments are stored, in job_args.yml

- pin_id:

  the unique identifier for the job

- event_type:

  either 'start' or 'end', returns the corresponding email for the given
  type

- job_results_url:

  the URL where the user can check the status of their job

## Value

the result of the sendmailR::sendmail() call
