# functions to measure:
# - process runtimes
# - handle log files

#' Stack rows of job logfiles
#' @param dir_job_results the base directory containing all job_results folders as subdirs
#' @return list of log table and a character vector of logfile paths that failed to read
# example
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   logs <- aggregate_logs(path_dev_results, verbose = TRUE)
aggregate_logs <- function(
  dir_job_results,
  latest_date = (Sys.Date() - 60), # default 60 days prior
  verbose = FALSE
) {

  # job results dirs
  vec_dir_results <- list.dirs(
    dir_job_results,
    recursive = FALSE,
    full.names = TRUE
  )
  # include only folders from `latest_date`
  vec_dir_results <- vapply(
    vec_dir_results,
    FUN = function(x) {
      lgl <- (as.integer(file.mtime(x))) < (as.integer(as.POSIXct(latest_date)))
      ifelse(lgl, yes = x, no = "")
    },
    FUN.VALUE = character(1)
  )
  vec_dir_results <- vec_dir_results[vec_dir_results != ""]

  # if verbose show warnings as they occur
  # and restore previous setting at end of function
  opt_warn <- getOption("warn")
  if (verbose) {options(warn = 1)}

  # empty df
  df_log <- tibble::tibble()
  # logfiles that could not be read 
  failed_reads <- c()
  # empty logsfiles
  empty_logs <- c()
  # log does not exist
  dne_logs <- c()
  for (dir in vec_dir_results) {
    tryCatch(
      expr = {
        path_logfile <- file.path(dir, "logfile.tsv")
        if (!file.exists(path_logfile)) {
          stringr::str_glue("log does not exist: {file.path(dir, 'logfile.tsv')}") |>
            warning()
          dne_logs <- append(x = dne_logs, values = dir)
          next
        }
        df <- readr::read_tsv(
          path_logfile,
          show_col_types = F
        )
        if (nrow(df) == 0L) {
          stringr::str_glue("empty log: {file.path(dir, 'logfile.tsv')}") |>
            warning()
          empty_logs <- append(x = empty_logs, values = dir)
          next
        }
        df_log <- df_log |> dplyr::bind_rows(df)
      },
      error = function(e) {
        failed_reads <- append(x = failed_reads, values = dir)
        if (verbose) {
          stringr::str_glue("failed to read log: {file.path(dir, 'logfile.tsv')}") |>
            warning()
        }
      }
    )
  }
  # reset warning option
  if (verbose) {options(warn = opt_warn)}
  return(
    list(
      "df_log" = df_log,
      "failed_reads" = failed_reads,
      "empty_logs" = empty_logs,
      "dne_logs" = dne_logs
    )
  )
}

#' calc statistic for processes in MolEvolvR logfiles
#' @param df_log `df_log` element from the return list of `aggregate_logs()`
#' @columns a single or multiple column used to calculate the statistic
#' @columns_group_by optionally group the data before statistic calculation
#' @metric a function to apply on the column(s)
#' @return if group_by is empty: a single value; else: tibble from summarise
# example
# median (no row grouping)
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   logs <- aggregate_logs(path_dev_results)
#   dblast_median <- logs$df_log |> calc_log_process_stat(columns = 'dblast', f = median)
# median (group_by query)
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   logs <- aggregate_logs(path_dev_results)
#   dblast_median_by_query <- logs$df_log |> calc_log_process_stat(columns = 'dblast', f = median, columns_group_by = 'query')
calc_log_process_stat <- function(
  df_log, columns, f, ..., columns_group_by = character(0)
) {
  # do not group rows
  if (purrr::is_empty(columns_group_by)) {
    result <- df_log |>
      tidyr::drop_na(columns) |>
      dplyr::mutate('stat' = dplyr::across(columns,  f, ...)) |>
      dplyr::pull('stat') |>
      dplyr::pull() |>
      unique()
  } else {
    # group rows
    result <- df_log |>
      tidyr::drop_na(columns) |>
      dplyr::group_by(.data[[columns_group_by]]) |>
      dplyr::mutate('stat' = dplyr::across(columns, f, ...)) |>
      dplyr::ungroup()
  }
  return(result)
}
