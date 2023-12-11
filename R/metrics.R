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

# use modification time of status.txt to estimate submission times
get_df_t_submit <- function(dir_job_results) {

  # job results dirs
  vec_dir_results <- list.dirs(
    dir_job_results,
    recursive = FALSE,
    full.names = TRUE
  )

  # search job folders for `status.txt`
  # use modification time to estimate submission date
  col_month <- character()
  col_year <- character()
  for (directory in vec_dir_results) {
    # catch unforeseen errors
    tryCatch(
      expr = {
        file_status_txt <- file.path(directory, "status.txt")
          if (file.exists(file_status_txt)) {
            modification_date <- file.mtime(file_status_txt)
            col_month <- append(
              x = col_month,
              values = modification_date |> months()
            )
            col_year <- append(
              x = col_year,
              values = modification_date |> format("%Y")
            )
          }
        },
      error = function(e) {
        msg <- stringr::str_glue("failed to measure submission date for '{directory}'")
        warning(msg)
      }
    )
  }
  df_t_submit <- tibble::tibble("month" = col_month, "year" = col_year)
  return(df_t_submit)

}




plot_df_t_submit <- function(df_t_submit) {
  df_n_submissions <- df_t_submit |> dplyr::group_by(month, year) |>
    dplyr::summarise(submissions = dplyr::n(), .groups = 'drop') |>
    dplyr::arrange(year, month)
  p <- ggplot2::ggplot(data = df_n_submissions) +
    ggplot2::aes(x = paste(year, month, sep = '-'), y = submissions) +
    ggplot2::theme_minimal() +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Date", y = "Number of submissions") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 30, face = 'bold'),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.text = ggplot2::element_text(size = 20),
      axis.title = ggplot2::element_text(size = 22, face = "bold"),
    )
  return(p)
}
path_dev_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
path_prod_results <- "/data/molevolvr_transfer/hpc-cluster-tests/job_results"
df_t_submit_dev <- get_df_t_submit(path_dev_results)
df_t_submit_prod <- get_df_t_submit(path_prod_results)
plot_df_t_submit(df_t_submit = df_t_submit_prod)