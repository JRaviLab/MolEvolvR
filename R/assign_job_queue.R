
#' Construct list where names (MolEvolvR input opts) point to processes
#'
#'
#' @return list where names (MolEvolvR input opts) point to processes
#'
#' example: list_opts2procs <- make_opts2procs
make_opts2procs <- function() {
  opts2processes <- list(
    "homology_search" = c("dblast", "dblast_cleanup"),
    "domain_architecture" = c("iprscan", "ipr2lineage", "ipr2da"),
    "any" = c("blast_clust", "clust2table") # processes always present agnostic of input opts 
  )
  return(opts2processes)
}

#' Use MolEvolvR input options to get associated processes
#'
#' @param input_opts character vector of MolEvolvR input opts
#'
#' @return character vector of process names that will execute given
#' the input options
#' 
#' example:
#' input_opts <- c("homology_search", "domain_architecture")
#' procs <- map_input_opts2procs(input_opts)
map_input_opts2procs <- function(input_opts) {
  # append 'any' to add procs that always run
  input_opts <- c(input_opts, "any")
  opts2proc <- make_opts2procs()
  # setup index for opts2proc based on input
  idx <- which(names(opts2proc) %in% input_opts)
  # extract processes that will run
  procs <- opts2proc[idx] |> unlist()
  return(procs)
}

#' Scrape logs and calculate median processes
#'
#' for now execute this code outside of container
#'
#' @return list: names are processes and point to a single value of 
#' median seconds for a given process
#'
#' see molevol_scripts/R/metrics.R for info on functions called here
#'
#' example:
#' input_opts <- c("homology_search", "domain_architecture")
#' procs <- map_input_opts2procs(input_opts)
get_proc_medians <- function(job_results_folder) {
  # first, use rprojroot to identify the molevol_scripts base folder so we can
  # do project-local imports
  common_root <- rprojroot::has_file(".molevol_root")

  source(common_root$find_file("molevol_scripts", "R", "metrics.R"))
  
  # aggregate logs from
  path_prod_results <- job_results_folder
  path_log_data <- common_root$find_file("molevol_scripts", "log_data", "prod_logs.rda")

  # ensure the folder exists to the location
  dir.create(dirname(path_log_data), recursive=TRUE)

  if (!file.exists(path_log_data)) {
    logs <- aggregate_logs(path_prod_results, latest_date = Sys.Date() - 60)
    save(logs, file = path_log_data)
  } else {
    load(path_log_data) # loads the logs object
  }
  df_log <- logs$df_log

  dblast_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'dblast', f = median)
  dblast_cleanup_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'dblast_cleanup', f = median)
  iprscan_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'iprscan', f = median)
  ipr2lineage_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'ipr2lineage', f = median)
  ipr2da_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'ipr2da', f = median)
  blast_clust_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'blast_clust', f = median)
  clust2table_median <- df_log |>
    calc_log_process_stat(na.rm = T, columns = 'clust2table', f = median)
  return(
    list(
      "dblast" = dblast_median,
      "dblast_cleanup" = dblast_cleanup_median,
      "iprscan" = iprscan_median,
      "ipr2lineage" = ipr2lineage_median,
      "ipr2da" = ipr2da_median,
      "blast_clust" = blast_clust_median,
      "clust2table" = clust2table_median
    )
  )
}

#' Write a table of 2 columns: 1) process and 2) median seconds
#'
#' @param filepath path to save table
#'
#' @return tibble: 2 columns: 1) process and 2) median seconds
#'
#' example: write_proc_medians_table(
#'   "/data/scratch/janani/molevolvr_out/",
#'   "/data/scratch/janani/molevolvr_out/log_tbl.tsv"
#' )
write_proc_medians_table <- function(
    job_results_folder,
    filepath
) {
proc_medians <- get_proc_medians(job_results_folder)
df_proc_medians <- proc_medians |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to = "process",
    values_to = "median_seconds"
  ) |>
  dplyr::arrange(dplyr::desc(median_seconds))
  readr::write_tsv(df_proc_medians, file = filepath)
  return(df_proc_medians)
}

#' Compute median process runtimes, then write a YAML list of the processes and
#' their median runtimes in seconds to the path specified by 'filepath'.
#' 
#' The default value of filepath is the value of the env var
#' MOLEVOLVR_PROC_WEIGHTS, which get_proc_weights() also uses as its default
#' read location.
#'
#' @param filepath path to save YAML file
#'
#' example: write_proc_medians_yml(
#'   "/data/scratch/janani/molevolvr_out/",
#'   "/data/scratch/janani/molevolvr_out/log_tbl.tsv"
#' )
write_proc_medians_yml <- function(
  job_results_folder,
  filepath=Sys.getenv("MOLEVOLVR_PROC_WEIGHTS", "/data/scratch/janani/molevolvr_out/job_proc_weights.yml")
) {
  medians <- get_proc_medians(job_results_folder)
  yaml::write_yaml(medians, filepath)
}

# molevolvr backend process weight list based on the median walltimes
# for each process (proc)
# example: get_proc_weights()
get_proc_weights <- function(
  medians_yml_path=Sys.getenv("MOLEVOLVR_PROC_WEIGHTS", "/data/scratch/janani/molevolvr_out/job_proc_weights.yml")
) {
  proc_weights <- tryCatch(
    {
      # attempt to read the weights from the YAML file produced by
      # write_proc_medians_yml()
      if (stringr::str_trim(medians_yml_path) == "") {
        stop(
          stringr::str_glue("medians_yml_path is empty ({medians_yml_path}), returning default weights")
        )
      }

      proc_weights <- yaml::read_yaml(medians_yml_path)
    },
    error=function(cond) {
      proc_weights <- list(
        "dblast" = 2810,
        "iprscan" = 1016,
        "dblast_cleanup" = 79,
        "ipr2lineage" = 18,
        "ipr2da" = 12,
        "blast_clust" = 2,
        "clust2table" = 2
      )
    }
  )

  return(proc_weights)
}

#' Given MolEvolvR input options and number of inputs,
#' calculate the total estimated walltime for the job
#'
#' @param input_opts character vector of MolEvolvR input options
#' (see make_opts2procs for the options)
#' @param n_inputs total number of input proteins
#'
#' @return total estimated number of seconds a job will process (walltime)
#'
#' example: input_opts2est_walltime(c("homology_search", "domain_architecture"), n_inputs = 3, n_hits = 50L)
input_opts2est_walltime <- function(input_opts, n_inputs = 1L, n_hits = NULL) {
  # to calculate est walltime for a homology search job, the number of hits
  # must be provided
  validation_fail <- is.null(n_hits) && "homology_search" %in% input_opts
  stopifnot(!validation_fail)

  proc_weights <- get_proc_weights()
  # sort process weights by names and convert to vec
  proc_weights <- proc_weights[order(names(proc_weights))] |> unlist()
  all_procs <- names(proc_weights) |> sort()
  # get processes from input options and sort by names
  procs_from_opts <- map_input_opts2procs(input_opts)
  procs_from_opts <- sort(procs_from_opts)
  # binary encode: yes proc will run (1); else 0
  binary_proc_vec  <- dplyr::if_else(all_procs %in% procs_from_opts, 1L, 0L)
  # dot product of weights and procs to run; scaled by the number of inputs
  est_walltime_queries <- (n_inputs * (binary_proc_vec %*% proc_weights)) |>
    as.numeric()
  # calculate the additional processes to run for the homologous hits
  if ("homology_search" %in% input_opts) {
    opts2procs <- make_opts2procs()
    # exclude the homology search processes for the homologous hits
    procs2exclude_for_homologs <- opts2procs[["homology_search"]]
    procs_homologs <- procs_from_opts[!(procs_from_opts %in% procs2exclude_for_homologs)]
    binary_proc_vec_homolog  <- dplyr::if_else(all_procs %in% procs_homologs, 1L, 0L)
    # add the estimated walltime for processes run on the homologous hits
    est_walltime_final <- est_walltime_queries + (n_hits * (binary_proc_vec_homolog %*% proc_weights))
  }
  return(est_walltime_final)
}

#' Decision function to assign job queue
#'
#' @param t_sec_estimate estimated number of seconds a job will process
#' (from `input_opts2est_walltime`)
#' @param t_long threshold value that defines the lower bound for assigning a
#' job to the "long queue"
#'
#' @return a string of "short" or "long"
#'
#' example:
#' input_opts2est_walltime(c("homology_search", "domain_architecture"), 3) |>
#'   assign_job_queue()
assign_job_queue <- function(
  t_sec_estimate,
  t_cutoff = 21600 # 6 hours
) {
  queue <- ifelse(t_sec_estimate > t_cutoff, "long", "short")
  return(queue)
}

#' Plot the estimated runtimes for different input options and number
#' of inputs
#'
#' this function was just for fun; very, very messy code
#'
#' @return line plot object
#'
#' example:
#' p <- plot_estimated_walltimes()
#' ggplot2::ggsave(filename = "/data/molevolvr_transfer/molevolvr_dev/molevol_scripts/docs/estimate_walltimes.png", plot = p)
plot_estimated_walltimes <- function() {
  opts <- make_opts2procs() |> names()
  # get all possible submission permutations (powerset)
  get_powerset <- function(vec) {
    # generate powerset (do not include empty set)
    n <- length(vec)
    indices <- 1:n
    powerset <- lapply(1:n, function(x) combn(indices, x, simplify = FALSE))
    powerset <- unlist(powerset, recursive = FALSE)
    powerset <- lapply(powerset, function(index) vec[index])
    powerset
  }
  opts_power_set <- get_powerset(opts)
  est_walltimes <- list()
  for (n_inputs in 1:20) {
    est_walltimes <- append(
      x = est_walltimes,
      values = sapply(
        opts_power_set,
        FUN = function(input_opts) {
          est_walltime <- input_opts2est_walltime(input_opts, n_inputs = n_inputs)
          names(est_walltime) <- paste0(input_opts, collapse = "_")
          est_walltime
        }
      ) 
    )
  }
  # concat all results to their unique names
  est_walltimes <- tapply(
    unlist(
      est_walltimes, use.names = FALSE),
      rep(names(est_walltimes),
      lengths(est_walltimes)
    ),
    FUN = c
  )
  df_walltimes <- est_walltimes |>
    unlist() |> 
    matrix(nrow = length(est_walltimes[[1]]), ncol = length(names(est_walltimes)))
  colnames(df_walltimes) <- names(est_walltimes)
  df_walltimes <- df_walltimes |> tibble::as_tibble()
  # rm any col or powerset outcome without the "any" processes
  col_idx_keep <- grep(pattern = "any$", x = names(df_walltimes))
  df_walltimes <- df_walltimes |>
    dplyr::select(col_idx_keep)
  # bind n_inputs
  df_walltimes <- df_walltimes |>
    dplyr::mutate(n_inputs = 1:20)
  df_walltimes <- tidyr::gather(df_walltimes, key = "input_opts", value = "est_walltime", -n_inputs)
  # sec to min
  df_walltimes <- df_walltimes |>
    dplyr::mutate(est_walltime = est_walltime / 60)
  p <- ggplot2::ggplot(df_walltimes, ggplot2::aes(x = n_inputs, y = est_walltime, color = input_opts)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "MolEvolvR estimated runtimes",
        x = "Number of Inputs",
        y = "Estimated walltime (log_10(minutes))")
  return(p)
}
