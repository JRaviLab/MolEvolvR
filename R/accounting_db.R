# library(RMariaDB)
# library(DBI)
# library(dplyr)

#' getCon
#'
#' @description
#' Gets a connection to the accounting db. You should
#' close this connection later via dbDisconnect(con).
#'
#' @importFrom DBI dbConnect
#' @importFrom RMariaDB MariaDB
#'
#' @return A DBI connection object
#' @export
#'
getCon <- function() {
    # get authentication info from environment variables
    MARIADB_USER <- Sys.getenv("MARIADB_USER")
    MARIADB_DATABASE <- Sys.getenv("MARIADB_DATABASE")
    MARIADB_ROOT_PASSWORD <- Sys.getenv("MARIADB_ROOT_PASSWORD")
    MARIADB_PASSWORD <- Sys.getenv("MARIADB_PASSWORD")

    con <- dbConnect(
        RMariaDB::MariaDB(),
        user = MARIADB_USER,
        password = MARIADB_PASSWORD,
        dbname = MARIADB_DATABASE,
        host = "accounting"
    )
}

#' getJobsForCode
#'
#' @description
#' Given a job code, retrieves SLURM jobs where that code occurs in their name
#' (i.e., the "job_name" column). Optionally, you may specify a vector of column
#' names to return; if unspecified, returns all columns.
#' 
#' @param code A job code to search for in job names
#'
#' @importFrom DBI dbGetQuery dbDisconnect dbQuoteIdentifier
#'
#' @return Data frame with 'cols' columns (or all columns if 'cols') of jobs
#' @export
#'
getJobsForCode <- function(code, cols=NA) {
    con <- getCon()

    colSpec <- if (!any(is.na(cols))) {
        toString(dbQuoteIdentifier(con, cols))
    } else { "*" }

    query <- paste0(
        "SELECT ", colSpec, " FROM localcluster_job_table WHERE job_name like ?"
    )

    result <- dbGetQuery(
        con, query,
        params=list(paste0("%", code, "%"))
    )

    dbDisconnect(con)

    return(result)
}

#' Get jobs with a specific SLURM job state.
#' 
#' States are stored as integers, whose values derive from the file
#' https://github.com/SchedMD/slurm/tree/master/slurm/slurm.h. The possible
#' values and enumeration values are as follows.
#' \itemize{
#'   \item 0:  JOB_PENDING: queued waiting for initiation
#'   \item 1:  JOB_RUNNING: allocated resources and executing
#'   \item 2:  JOB_SUSPENDED: allocated resources, execution suspended
#'   \item 3:  JOB_COMPLETE: completed execution successfully
#'   \item 4:  JOB_CANCELLED: cancelled by user
#'   \item 5:  JOB_FAILED: completed execution unsuccessfully
#'   \item 6:  JOB_TIMEOUT: terminated on reaching time limit
#'   \item 7:  JOB_NODE_FAIL: terminated on node failure
#'   \item 8:  JOB_PREEMPTED: terminated due to preemption
#'   \item 9:  JOB_BOOT_FAIL: terminated due to node boot failure
#'   \item 10: JOB_DEADLINE: terminated on deadline
#'   \item 11: JOB_OOM: experienced out of memory error
#' }
#' (See https://slurm.schedmd.com/squeue.html#SECTION_JOB-STATE-CODES for
#' the meanings of these codes.)
#' 
#' For example, to search for all running, suspended, and pending jobs, you'd
#' provide the states argument states=c(0, 1, 2). If you wanted just running
#' jobs, you'd specify states=1.
#' 
#' @param cols Column names to return; if unspecified, returns all columns.
#' @param states integers for states to search for; you can specify either a single
#' integer or a vector. if unspecified, returns all rows.
#'
#' @importFrom DBI dbDisconnect
#' @importFrom dplyr any_of collect filter tbl select
#'
#' @return Data frame with 'cols' columns (or all columns if 'cols') of jobs
#' with the specified states, or all jobs if 'states' is unspecified
#' @export
#'
getJobsWithStates <- function(states=NA, cols=NA) {
    con <- getCon()

    # produces a "lazy" query object, i.e. the description
    # of a query; it's not materialized into an actual value
    # until it's used, or explicitly materialized by collect(<lazy_query>).
    # since it's just a description, it can be constrained later on
    # using %>% with additional filter(), select(), etc. calls.
    lazy_result <- tbl(con, "localcluster_job_table")

    # if states were specified, filter the jobs by the states
    # (leaving the select() out implies we want all jobs)
    if (!any(is.na(states))) {
        lazy_result <- lazy_result %>% filter(state %in% states)
    }

    # if columns were specified, apply those via select()
    # (leaving the select() out implies we want all columns)
    if (!any(is.na(cols))) {
        lazy_result <- lazy_result %>% select(any_of(cols))
    }

    # executes the query, producing an actual dataframes
    # (this must be done before the connection's closed)
    result_df <- collect(lazy_result)

    dbDisconnect(con)

    return(result_df)
}
