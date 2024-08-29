# Description
# functions for sending emails to users for various stages of job processing.
# currently includes job start and job completion emails.
# uses sendmailR (https://cran.r-project.org/web/packages/sendmailR/index.html)
# under the hood to actually send the emails
# Usage
# 1. source("job_status_emails.R")
# 2. call:
#     # event_type can be 'start' or 'end'
#     send_job_status_email(notify_email, job_dir, pin_id, event_type)
# Return
#   unfortunately, there is no return value for the underlying sendmailR methods
#
# library(sendmailR, warn.conflicts=F, quietly=T)
# library(shiny, warn.conflicts=F, quietly=T) # includes tags that we use for HTML-formatting message

#' Given a pin_id, returns the URL where the user can check the status of their job
#'
#' The 'base_url' parameter is set from the environment variable 'BASE_URL' if it's
#' available, and defaults to "http://jravilab.org/molevolvr/" if it's not.
#'
#' @param pin_id the unique identifier for the job
#' @param base_url the base URL for the MolEvolvR web app
#'
#' @return the URL where the user can check the status of their job
#' @export
#'
make_job_results_url <- function(
        pin_id,
        base_url = Sys.getenv("BASE_URL", unset = "http://jravilab.org/molevolvr/")) {
    return(paste0(base_url, "?r=", pin_id, "&p=home"))
}

# stores human-readable names for fields and their possible values
# (if an identifier isn't found here, it should be used verbatim)
# fields_metadata <- list(
#     submission_type=list(
#         name="Submission Type",
#         values=list(
#             'full'="FASTA",
#             'phylo'="Phylogenetic Analysis",
#             'da'="Domain Architecture",
#             'blast'="BLAST",
#             'dblast'="BLAST"
#         )
#     ),
#     homology_search=list(
#         name="Homology Search?"
#     ),
#     submitter_email=list(
#         name="Submitter Email"
#     ),
#     database=list(
#         name="Database for homology search"
#     ),
#     nhits=list(
#         name="Maximum hits"
#     ),
#     evalue=list(
#         name="E-value cutoff"
#     ),
#     includes_ncbi_acc=list(
#         name="Contains NCBI accessions?",
#         values=list(
#             "TRUE"="Yes",
#             "FALSE"="No"
#         )
#     ),
#     advanced_options=list(
#         name="Advanced Options",
#         values=list(
#             "domain_architecture"="Domain Architecture",
#             "homology_search"="Homology Search",
#             "phylogenetic_analysis"="Phylogenetic Analysis"
#         )
#     ),
#     job_code=list(
#         name="Job Code"
#     )
# )

#' Format job arguments into html-formatted key/value pairs, for including
#' in an email
#'
#' @param job_args
#' a list of job arguments, e.g. as read from the job_args.yml file
#'
#' @return
#' a list of HTML-formatted key/value pairs
#' @export
#'
#' @examples
#' \dontrun{
#' format_job_args("/data/scratch/janani/molevolvr_out/Ba5sV1_full")
#' }
format_job_args <- function(job_args) {
    # format job arguments into html-formatted key/value pairs
    job_args_list <- tags$ul(lapply(names(job_args), function(key) {
        # look up human labels for field names, values, if available
        # (if not, just use the keys and values as-is)
        field_meta <- fields_metadata[[key]]
        human_key <- tryCatch(
            {
                toString(
                    ifelse(!is.null(field_meta), field_meta$name, key)
                )
            },
            error = function(e) {
                return(value)
            }
        )

        value <- job_args[[key]]

        # this function takes a value and uses the field metadata to
        # look up a human-readable version of the value, if available.
        # otherwise, it returns the original vlaue.
        humanize_value <- function(value) {
            tryCatch(
                {
                    ifelse(
                        !is.null(field_meta) && !is.null(field_meta$values) && !is.null(field_meta$values[[toString(value)]]),
                        field_meta$values[[value]], value
                    )
                },
                error = function(e) {
                    return(value)
                }
            )
        }

        # invoking lappy(humanize) deals with both single values and lists of
        # values, which is necessary because, e.g., the advanced_options field
        # is a list.
        # the result of running toString() on a list is a comma-delimited list,
        # which works fine for our purposes.
        human_value <- toString(lapply(value, humanize_value))

        tags$li(
            tags$b(paste0(human_key, ":")),
            human_value
        )
    }))

    return(job_args_list)
}

#' Produces a mail message that can be sent to a user when their job is accepted.
#' Used by the send_job_status_email() method.
#'
#' @param job_dir
#' the directory where the job's arguments are stored, in job_args.yml
#' @param pin_id
#' the unique identifier for the job
#' @param job_results_url
#' the URL where the user can check the status of their job
#' @param event_type
#' either 'start' or 'end', returns the corresponding email for the given type
#' @param context
#' a list of additional values, e.g. job runtime info, that can be used in the template emails
#'
#' @importFrom htmltools htmlTemplate
#' @importFrom sendmailR mime_part
#' @importFrom yaml read_yaml
#'
#' @return
#' the result of the sendmailR::sendmail() call
#' @export
get_job_message <- function(job_dir, pin_id, job_results_url, event_type, context) {
    # pull the set of args written to dir/job_args.yml, so we
    # can send it in the email
    job_args <- yaml::read_yaml(file.path(job_dir, "job_args.yml"))
    job_args_list <- format_job_args(job_args)

    # determine which template to use based on the event type
    if (event_type == "start") {
        template <- "/data/research/jravilab/molevol_scripts/templates/job_start_email/job_start_email.html"
    } else if (event_type == "end") {
        template <- "/data/research/jravilab/molevol_scripts/templates/job_end_email/job_end_email.html"
    } else {
        stop("Invalid event type (expected 'start' or 'end'): ", event_type)
    }

    # use the vars from above injected into a template to generate the final email text
    message <- mime_part(
        paste0(
            htmlTemplate(
                template,
                pin_id = pin_id,
                job_results_url = job_results_url,
                job_args_list = job_args_list,
                context = context
            )
        ),
        type = "text/html"
    )

    return(message)
}

#' Sends a "job accepted" email to a user when their job is accepted,
#' including details about the job submission and how to check its status.
#'
#' @param notify_email
#' the email address to send the notification to
#' @param job_dir
#' the directory where the job's arguments are stored, in job_args.yml
#' @param pin_id
#' the unique identifier for the job
#' @param job_results_url
#' the URL where the user can check the status of their job
#' @param event_type
#' either 'start' or 'end', returns the corresponding email for the given type
#'
#' @importFrom sendmailR sendmail
#'
#' @return
#' the result of the sendmailR::sendmail() call
#' @export
send_job_status_email <- function(notify_email, job_dir, pin_id, event_type, context = NULL) {
    # -------------------------------------------------
    # --- step 1. build the email subject and contents
    # -------------------------------------------------

    # determine the subject based on the event type
    if (event_type == "start") {
        mail_subject <- paste0("MolEvolvR Job Accepted (ID: ", pin_id, ")")
    } else if (event_type == "end") {
        mail_subject <- paste0("MolEvolvR Job Completed (ID: ", pin_id, ")")
    } else {
        stop("Invalid event type (expected 'start' or 'end'): ", event_type)
    }

    # construct the job results URL from the pin_id
    job_results_url <- make_job_results_url(pin_id)

    # produce a formatted email message from the arguments and template
    message <- get_job_message(
        job_dir, pin_id, job_results_url, event_type, context
    )

    # -------------------------------------------------
    # --- step 2. send the email via sendmailr or curl
    # -------------------------------------------------

    # if verbose, output emails
    mail_verbose <- (Sys.getenv("SMTP_VERBOSE", unset = 1) == 1)

    if (Sys.getenv("SMTP_USE_CURL", unset = 1) == 1) {
        # when using curl, we need a full URL to the server of the form
        # smtp://<host>:<port>

        # note that if the port is 25, for some reason sendmailR will append
        # it again and we'll end up with :25:25, which is invalid. instead, we
        # have to check that it's not 25 and then append the port...

        target_port <- Sys.getenv("SMTP_SSL_PORT")

        if (target_port != "25") {
            target_server <- paste0("smtp://", Sys.getenv("SMTP_HOST"), ":", target_port)
        } else {
            target_server <- paste0("smtp://", Sys.getenv("SMTP_HOST"))
        }

        # and finally send the email using cURL
        result <- sendmail(
            from = Sys.getenv("SMTP_EMAIL"),
            to = c(notify_email),
            subject = mail_subject,
            msg = message,
            engine = "curl",
            engineopts = list(
                username = Sys.getenv("SMTP_EMAIL"),
                password = Sys.getenv("SMTP_PASSWORD")
            ),
            control = list(
                smtpServer = target_server,
                verbose = mail_verbose
            )
        )
    } else {
        # when not using curl, we don't need a URL, just the host
        # note that CU's SMTP server balks at STARTTLS, so we have to use
        # the non-curl implementation
        result <- sendmail(
            from = Sys.getenv("SMTP_EMAIL"),
            to = c(notify_email),
            subject = mail_subject,
            msg = message,
            control = list(
                smtpServer = Sys.getenv("SMTP_HOST"),
                verbose = mail_verbose
            )
        )
    }
    return(result)
}
