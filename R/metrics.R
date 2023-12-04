# functions to measure:
# - process runtimes
# - handle log files

#' Stack rows of job logfiles
#' @param dir_job_results the base directory containing all job_results folders as subdirs
# example
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   df_logs <- aggregate_logs(path_dev_results, verbose = TRUE)
aggregate_logs <- function(dir_job_results, verbose = FALSE) {
  vec_dir_results <- list.dirs(
    dir_job_results,
    recursive = FALSE,
    full.names = TRUE
  )
  # if verbose show warnings as they occur
  # and restore previous setting at end of function
  opt_warn <- getOption("warn")
  if (verbose) {options(warn = 1)}

  # empty df
  df_log <- tibble::tibble()
  # logfiles that could not be read (or they were never written)
  failed_reads <- c()
  for (dir in vec_dir_results) {
    tryCatch(
      expr = {
        df <- readr::read_tsv(
          file.path(dir, "logfile.tsv"),
          show_col_types = F
        )
        df_log <- df_log |> dplyr::bind_rows(df)
      },
      error = function(e) {
        failed_reads <- append(
          failed_reads,
          dir
        )
        if (verbose) {
          stringr::str_glue("failed to read {file.path(dir, 'logfile.tsv')}") |>
            warning()
        }
      }
    )
  }
  # reset warning option
  if (verbose) {options(warn = opt_warn)}
  return(list("df_log" = df_log, "failed_reads" = failed_reads))
}

#' calc statistic for processes in MolEvolvR logfiles
#' @param df_log `df_log` element from the return list of `aggregate_logs()`
#' @columns a single or multiple column used to calculate the statistic
#' @columns_group_by optionally group the data before statistic calculation
#' @metric a function to apply on the column(s)
#' @return if group_by is empty: a single value; else: tibble from summarise
# example
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   df_log <- aggregate_logs(path_dev_results)
#   result <- df_log$df_log |> calc_log_process_stat("dblast", f = mean)
# example group_by
#   path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#   df_log <- aggregate_logs(path_dev_results)
#   df_log$df_log |> count_log_process_stat(columns = 'dblast', columns_group_by = 'query')
calc_log_process_stat <- function(
  df_log, columns, f, ..., columns_group_by = character(0)
) {
  # do not group rows
  if (purrr::is_empty(columns_group_by)) {
    result <- df_log |>
      tidyr::drop_na(columns) |>
      dplyr::mutate('stat' = dplyr::across(columns,  f)) |>
      dplyr::pull('stat') |>
      dplyr::pull() |>
      unique()
  } else {
    # group rows
    result <- df_log |>
      tidyr::drop_na(columns) |>
      dplyr::group_by(columns_group_by) |>
      dplyr::summarise('stat' = f)
  }
  return(result)
}

