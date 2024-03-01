#!/usr/bin/env Rscript

# This script sends an email to the submitter of a job to notify them that the
# job has completed. It is intended to be called by 00_submit_full.R and should
# not be called directly.

# Arguments:
# 1. job_dir: the directory where the job results are stored; must contain
#    job_args.yml
# 2. submitter_email: the email address to send the notification to
# 3. pin_id: the unique identifier for the job
# 4. dep_jobs: a comma-delimited list of job IDs that this job depends on
#    (these jobs are the ones that actually did the work; we query the
#    accounting subsystem to get more information about their execution and
#    resource usage)

# we need the send_job_status_email() method from this file
source("/data/research/jravilab/molevol_scripts/R/job_status_emails.R")

# ------------------------------------------------------------------------
# --- supporting methods
# ------------------------------------------------------------------------

# true if v is null or an empty string
null_or_empty <- function(v) { is.null(v) || stringr::str_trim(v) == "" }

# value if not null/NA, else unset (default '')
val_or_else <- function(v, unset='') { if (is.null(v) || is.na(v)) { unset } else { v } }

# converts a number with a size specifier, e.g. 32000M, to bytes
parse_size <- function(size) {
    if (is.null(size) || size == "") return(NA)

    if (grep("K$", size)) {
        return(as.integer(sub("K$", "", size)) * 1024)
    } else if (grepl("M$", size)) {
        return(as.integer(sub("M$", "", size)) * 1024^2)
    } else if (grepl("G$", size)) {
        return(as.integer(sub("G$", "", size)) * 1024^3)
    } else if (grepl("T$", size)) {
        return(as.integer(sub("T$", "", size)) * 1024^4)
    } else {
        # just get the number part if the size isn't recognized or present at all
        return(as.integer(stringr::str_extract(size, "[0-9]+")))
    }
}

# runs a command and, if verbose, prints the command and its output
run_cmd <- function(cmd, verbose=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE) {
    if (verbose) write(paste0("Running command: ", cmd), file=stderr())
    
    result <- system(
        cmd, intern = TRUE,
        ignore.stdout = ignore.stdout,
        ignore.stderr = ignore.stderr
    )

    if (!is.null(attr(result, 'status'))) {
        write(paste0("Command failed with status ", attr(result, 'status')), file=stderr())
    }

    if (verbose) {
        write("Result: ", file=stderr())
        writeLines(result, con=stderr())
        write("", file=stderr())
    }

    return(result)
}

#' Get information about a set of finished jobs from the accounting system.
#' (Getting information about running jobs is possible via scontrol, but
#'  a bit too complicated to implement now, since we don't need it yet.)
#' 
#' @param job_ids an integer vector of job IDs to query
#' @param verbose if TRUE, print the job info to stderr
#' 
#' @return a list of the form: list(
#'  <job_id>=list(<job_field>=<value>, <...more fields>), <...more jobs>
#' )
#' 
get_completed_job_info <- function(job_ids, verbose=TRUE) {
    tryCatch({
        # store info about jobs by job_id in this list
        jobs <- list()

        # list the fields we want to pull for the job from sacct
        # (see https://slurm.schedmd.com/sacct.html#OPT_helpformat for options)
        fields <- c(
            "JobId", # required for filtering, excluded from email
            "JobName",
            "User",
            "Group",
            "Partition",
            "Start",
            "End",
            "State",
            "ReqMem",
            "MaxRSS",
            "NCPUS",
            "CPUTimeRaw",
            "TotalCPU",
            "NNodes",
            "WorkDir",
            "Submit",
            "SubmitLine",
            "Elapsed",
            "ExitCode",
            "AdminComment",
            "Comment",
            "Cluster",
            "NodeList",
            "TimeLimit",
            "TimelimitRaw",
            "JobIdRaw"
        )

        # loop over the jobs and query sacct for more information
        for (job_id in job_ids) {
            # get the job's accounting information via sacct
            cmd <- stringr::str_glue("sacct -j {job_id} --format={paste(fields, collapse=\",\")} --parsable2 --noheader")
            result <- run_cmd(cmd, verbose)

            # note that we'll get records back for each jobstep in the job,
            # and in the case of a batch job this is three records:
            # 1. the 'allocation' record in which the job is set up
            # 2. the 'batch' record where computation actually occurs
            # 3. the 'extern' record which finalizes the job
            # most of the job submission info is in the allocation step,
            # but we also need the batch step's MaxRSS to get the memory usage
            # (see https://slurm.schedmd.com/job_launch.html#job_record for more)

            # sacct produces per jobstep one line of the form
            # "<field1>|<field2>|...|<fieldN>", so we split each on the pipe
            all_jobsteps <- stringr::str_split(
                stringr::str_trim(result), "\\|",
                simplify=TRUE
            )

            # take the allocation step, i.e. the first one, as the main source of info
            job_info <- all_jobsteps[1,]
            # annotate the job_info with the MaxRSS column from the batch step,
            # with a little formatting to make it easier to read
            max_rss_idx <- match("MaxRSS", fields)
            formatted_maxrss <- utils:::format.object_size(parse_size(all_jobsteps[2, max_rss_idx]), units = "auto")
            job_info[max_rss_idx] <- formatted_maxrss # all_jobsteps[2, max_rss_idx]

            if (verbose) { writeLines(yaml::as.yaml(job_info), con=stderr()) }

            # if the job is not found, just skip it
            if (length(job_info) == 1) {
                write(paste0("Couldn't find job ", job_id, " in sacc, continuing to next"), file=stderr())
                next
            }

            # zip field names and values returned from sacct into a list
            job_info_list <- as.list(setNames(job_info, fields))

            # add a human-readable version of the id
            job_info_list$JobIdHuman <- stringr::str_replace(
                job_info_list$JobId, "_[0-9]+$", " (Array)"
            )

            # add this to the jobs
            jobs[[toString(job_info_list$JobId)]] <- job_info_list
        }

        return(jobs)
    }, error=function(e) {
        print(e)
        return(NULL)
    })
}

# ------------------------------------------------------------------------
# --- entrypoint
# ------------------------------------------------------------------------

# if not interactive (i.e. sourced from an R session), treat it like a shell script
if (!interactive()) {
    # unpack arguments to script
    args <- commandArgs(trailingOnly = TRUE)
    job_dir <- val_or_else(args[1]) # the folder in which job results are located
    submitter_email <- val_or_else(args[2]) # the email address of the submitter
    pin_id <- val_or_else(args[3]) # the job code, e.g. Ba5c8x
    dep_jobs <- val_or_else(args[4]) # a comma-delimited list of job IDs that this job depends on

    # if any of the arguments are empty, raise an error
    if (any(sapply(c(job_dir, submitter_email, pin_id), null_or_empty))) {
        write("Usage: 99_send_completion_mail.R <job_dir> <submitter_email> <pin_id> <dep_job_ids>", stderr())
        quit(status=1)
    }

    write(stringr::str_glue("* Sending completion email for job {pin_id} to {submitter_email}..."), file=stderr())

    # split dep_jobs on commas into vector of integers
    # and query for info about those jobs from the accounting system
    dep_job_ids_vec <- as.integer(unlist(strsplit(dep_jobs, ",")))
    total_job_info <- get_job_info(dep_job_ids_vec, verbose=FALSE)

    # get the minimum and maximum date values from total_job_info, i.e.
    # the minimum of the field "Start" and the maximum of the field "End"
    min_start <- lubridate::as_datetime(
        min(sapply(total_job_info, function(job) as.POSIXct(job$Start, format="%Y-%m-%dT%H:%M:%S")))
    )
    max_end   <- lubridate::as_datetime(
        max(sapply(total_job_info, function(job) as.POSIXct(job$End,   format="%Y-%m-%dT%H:%M:%S")))
    )

    # get a human-readable duration from min_start and max_end
    duration <- toString(hms::as_hms(difftime(max_end, min_start)))

    # convert min_start and max_end back into human-readable strings
    local_tz <- Sys.getenv("LOCAL_TIMEZONE", unset="America/Denver")
    outFormat <- "%Y-%m-%d %I:%M %p %Z"
    min_start_str <- format.POSIXct(min_start, outFormat, tz = local_tz)
    max_end_str   <- format.POSIXct(max_end,   outFormat, tz = local_tz)

    write(paste0("min_start: ", min_start_str), file=stderr())
    write(paste0("max_end: ", max_end_str), file=stderr())

    # currently we don't include information about the job that sends the email,
    # i.e., the currently running job in when this script is invoked, but if we
    # were to you'd get its ID here via Sys.getenv("SLURM_JOB_ID")

    send_job_status_email(
        submitter_email, job_dir, pin_id, 'end',
        context=list(
            duration=duration,
            min_start=min_start_str,
            max_end=max_end_str,
            jobs=total_job_info
        )
    )
}
