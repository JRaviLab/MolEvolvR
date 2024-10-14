# for now, we're using an env var, COMMON_SRC_ROOT, to specify this folder since
# the working directory is changed in many parts of the current molevolvr
# pipeline.
# to use this, construct paths like so: file.path(common_root, "path", "to", "file.R")
# for example, the reference for this file would be:
# file.path(common_root, "molevol_scripts", "R", "assignJobQueue.R")
common_root <- Sys.getenv("COMMON_SRC_ROOT")

#' Construct list where names (MolEvolvR advanced options) point to processes
#'
#' @return list where names (MolEvolvR advanced options) point to processes
#'
#' example: list_opts2procs <- mapOption2Process
#' @export
mapOption2Process <- function() {
  opts2processes <- list(
    "homology_search" = c("dblast", "dblast_cleanup"),
    "domain_architecture" = c("iprscan", "ipr2lineage", "ipr2da"),
    # processes always present agnostic of advanced options
    "always" = c("blast_clust", "clust2table")
  )
  return(opts2processes)
}

#' Use MolEvolvR advanced options to get associated processes
#'
#' @param advanced_opts character vector of MolEvolvR advanced options
#'
#' @return character vector of process names that will execute given
#' the advanced options
#'
#' example:
#' advanced_opts <- c("homology_search", "domain_architecture")
#' procs <- mapAdvOption2Process(advanced_opts)
#' @export
mapAdvOption2Process <- function(advanced_opts) {
  # append 'always' to add procs that always run
  advanced_opts <- c(advanced_opts, "always")
  opts2proc <- mapOption2Process()
  # setup index for opts2proc based on advanced options
  idx <- which(names(opts2proc) %in% advanced_opts)
  # extract processes that will run
  procs <- opts2proc[idx] |> unlist()
  return(procs)
}

#' Scrape MolEvolvR logs and calculate median processes
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' directory
#'
#' @importFrom dplyr across everything select summarise
#'
#' @return [list] names: processes; values: median runtime (seconds)
#'
#' see molevol_scripts/R/metrics.R for info on functions called here
#'
#' examples:
#'
#' 1)
#' dir_job_results <- "/data/scratch/janani/molevolvr_out"
#' list_proc_medians <- calculateProcessRuntime(dir_job_results)
#'
#' 2) from outside container environment
#' common_root <- "/data/molevolvr_transfer/molevolvr_dev"
#' dir_job_results <- "/data/molevolvr_transfer/molevolvr_dev/job_results"
#' list_proc_medians <- calculateProcessRuntime(dir_job_results)
#' @export
calculateProcessRuntime <- function(dir_job_results) {
  source(file.path(common_root, "molevol_scripts", "R", "metrics.R"))

  # aggregate logs from
  path_log_data <- file.path(common_root,
                              "molevol_scripts", "log_data", "prod_logs.rda")

  # ensure the folder exists to the location
  if (!dir.exists(path_log_data)) {
    dir.create(dirname(path_log_data),
                recursive = TRUE, showWarnings = FALSE)
  }

  # attempt to load pre-generated logdata
  if (!file.exists(path_log_data)) {
    logs <- aggregate_logs(dir_job_results, latest_date = Sys.Date() - 60)
    save(logs, file = path_log_data)
  } else {
    load(path_log_data) # loads the logs object
  }
  df_log <- logs$df_log
  procs <- c(
    "dblast", "dblast_cleanup", "iprscan",
    "ipr2lineage", "ipr2da", "blast_clust",
    "clust2table"
  )
  list_proc_medians <- df_log |>
    dplyr::select(dplyr::all_of(procs)) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        \(x) median(x, na.rm = TRUE)
      )
    ) |>
    as.list()
  return(list_proc_medians)
}

#' Write a table of 2 columns: 1) process and 2) median seconds
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' @param filepath path to save tsv file
#'
#' @importFrom dplyr arrange desc everything
#' @importFrom tibble as_tibble
#' @importFrom readr write_tsv
#' @importFrom tidyr pivot_longer
#'
#' @return [tbl_df] 2 columns: 1) process and 2) median seconds
#'
#' example: writeProcessRuntime2TSV(
#'   "/data/scratch/janani/molevolvr_out/",
#'   "/data/scratch/janani/molevolvr_out/log_tbl.tsv"
#' )
#' @export
writeProcessRuntime2TSV <- function(dir_job_results, filepath) {
  df_proc_medians <- calculateProcessRuntime(dir_job_results) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "process",
      values_to = "median_seconds"
    ) |>
    dplyr::arrange(dplyr::desc(median_seconds))

  # Write the resulting tibble to a TSV file
  readr::write_tsv(df_proc_medians, file = filepath)
  return(df_proc_medians)
}

#' Compute median process runtimes, then write a YAML list of the processes and
#' their median runtimes in seconds to the path specified by 'filepath'.
#'
#' The default value of filepath is the value of the env var
#' MOLEVOLVR_PROC_WEIGHTS, which writeProcessRuntime2YML() also uses as its default
#' read location.
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results directory
#' @param filepath [chr] path to save YAML file; if NULL, 
#'                 uses ./molevol_scripts/log_data/job_proc_weights.yml
#'
#' @importFrom yaml write_yaml
#'
#' @examples
#' \dontrun{
#' writeProcessRuntime2YML(
#'     "/data/scratch/janani/molevolvr_out/",
#'     "/data/scratch/janani/molevolvr_out/log_tbl.yml"
#' )
#' }
#' @export
writeProcessRuntime2YML <- function(dir_job_results, filepath = NULL) {
  medians <- calculateProcessRuntime(dir_job_results)
  yaml::write_yaml(medians, filepath)
}

#' Quickly get the runtime weights for MolEvolvR backend processes
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' directory
#'
#' @importFrom stringr str_glue str_trim
#' @importFrom yaml read_yaml
#'
#' @return [list] names: processes; values: median runtime (seconds)
#'
#' example: writeProcessRuntime2YML()
#' @export
getProcessRuntimeWeights <- function(medians_yml_path = NULL) {
  if (is.null(medians_yml_path)) {
    medians_yml_path <- file.path(common_root,
                                  "molevol_scripts",
                                  "log_data",
                                  "job_proc_weights.yml")
  }

  proc_weights <- tryCatch({
    # attempt to read the weights from the YAML file produced by
    # writeProcessRuntime2YML()
    if (stringr::str_trim(medians_yml_path) == "") {
      stop(
        stringr::str_glue("medians_yml_path is empty 
                          ({medians_yml_path}), returning default weights")
      )
    }

    proc_weights <- yaml::read_yaml(medians_yml_path)
  },
  # to avoid fatal errors in reading the proc weights yaml,
  # some median process runtimes have been hardcoded based on
  # the result of calculateProcessRuntime() from Jan 2024
  error = function(cond) {
    proc_weights <- list(
      "dblast" = 2810,
      "iprscan" = 1016,
      "dblast_cleanup" = 79,
      "ipr2lineage" = 18,
      "ipr2da" = 12,
      "blast_clust" = 2,
      "clust2table" = 2
    )
    proc_weights
  })

  return(proc_weights)
}

#' Given MolEvolvR advanced options and number of inputs,
#' calculate the total estimated walltime for the job
#'
#' @param advanced_opts character vector of MolEvolvR advanced options
#' (see mapOption2Process for the options)
#' @param n_inputs total number of input proteins
#'
#' @importFrom dplyr if_else
#' @importFrom stringr str_glue
#'
#' @return total estimated number of seconds a job will process (walltime)
#'
#' example: calculateEstimatedWallTimeFromOpts	(c("homology_search",
#'                                       "domain_architecture"),
#'                                       n_inputs = 3, n_hits = 50L)
#' @export
calculateEstimatedWallTimeFromOpts	 <- function(advanced_opts,
                                                  n_inputs = 1L,
                                                  n_hits = NULL,
                                                  verbose = FALSE) {
  # to calculate est walltime for a homology search job, the number of hits
  # must be provided
  validation_fail <- is.null(n_hits) && "homology_search" %in% advanced_opts
  stopifnot(!validation_fail)

  # Get process weights
  proc_weights <- writeProcessRuntime2YML()

  # sort process weights by names and convert to vec
  proc_weights <- proc_weights[order(names(proc_weights))] |> unlist()
  all_procs <- names(proc_weights) |> sort()
  # get processes from advanced options and sort by names
  procs_from_opts <- mapAdvOption2Process(advanced_opts)
  procs_from_opts <- sort(procs_from_opts)
  # binary encode: yes proc will run (1); else 0
  binary_proc_vec <- dplyr::if_else(all_procs %in% procs_from_opts, 1L, 0L)
  # dot product of weights and procs to run; scaled by the number of inputs
  est_walltime <- (n_inputs * (binary_proc_vec %*% proc_weights)) |>
    as.numeric()
  # calculate the additional processes to run for the homologous hits
  if ("homology_search" %in% advanced_opts) {
    opts2procs <- mapOption2Process()
    # exclude the homology search processes for the homologous hits
    procs2exclude_for_homologs <- opts2procs[["homology_search"]]
    procs_homologs <- procs_from_opts[!(procs_from_opts 
                                        %in% procs2exclude_for_homologs)]
    binary_proc_vec_homolog <- dplyr::if_else(all_procs 
                                              %in% procs_homologs, 1L, 0L)
    # add the estimated walltime for processes run on the homologous hits
    est_walltime <- est_walltime +
      (n_hits * (binary_proc_vec_homolog
                  %*% proc_weights) |> as.numeric())
  }
  if (verbose) {
    msg <- stringr::str_glue(
      "warnings from calculateEstimatedWallTimeFromOpts	():\n",
      "\tn_inputs={n_inputs}\n",
      "\tn_hits={ifelse(is.null(n_hits), 'null', n_hits)}\n",
      "\test_walltime={est_walltime}\n\n"
    )
    cat(file = stderr(), msg)
  }
  return(est_walltime)
}


#' Decision function to assign job queue
#'
#' @param t_sec_estimate estimated number of seconds a job will process
#' (from calculateEstimatedWallTimeFromOpts	())
#' @param t_long threshold value that defines the lower bound for assigning a
#' job to the "long queue"
#'
#' @return a string of "short" or "long"
#'
#' example:
#' calculateEstimatedWallTimeFromOpts	(c("homology_search",
#'                                         "domain_architecture"), 3) |>
#'   assignJobQueue()
#' @export
assignJobQueue <- function(
  t_sec_estimate,
  t_cutoff = 21600 # 6 hours
) {
  queue <- ifelse(t_sec_estimate > t_cutoff, "long", "short")
  return(queue)
}

#' Plot the estimated runtimes for different advanced options and number
#' of inputs
#'
#' this function was just for fun; very, very messy code
#'
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 aes geom_line ggplot labs
#' @importFrom tibble as_tibble
#'
#' @return line plot object
#'
#' example:
#' p <- plotEstimatedWallTimes()
#' ggplot2::ggsave(filename = "/data/molevolvr_transfer/molevolvr_
#'                 dev/molevol_scripts/docs/estimate_walltimes.png", plot = p)
#' @export
plotEstimatedWallTimes <- function() {
    opts <- mapOption2Process() |> names()
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
    for (i in 1:20) {
      est_walltimes <- append(
        x = est_walltimes,
        values = sapply(
          opts_power_set,
          FUN = function(advanced_opts) {
            # for simplicity, assume the default number of homologus hits (100)
            n_hits <- if ("homology_search" %in% advanced_opts) {
              100
            } else {
                NULL
              }
            est_walltime <- calculateEstimatedWallTimeFromOpts	(
              advanced_opts,
              n_inputs = i,
              n_hits = n_hits,
              verbose = TRUE
            )
            names(est_walltime) <- paste0(advanced_opts, collapse = "_")
            est_walltime
          }
        )
      )
    }
    # concat all results to their unique names
    est_walltimes <- tapply(
      unlist(
        est_walltimes,
        use.names = FALSE
      ),
      rep(
        names(est_walltimes),
        lengths(est_walltimes)
      ),
      FUN = c
    )
    df_walltimes <- est_walltimes |>
      unlist() |>
      matrix(nrow = length(est_walltimes[[1]]),
             ncol = length(names(est_walltimes)))
    colnames(df_walltimes) <- names(est_walltimes)
    df_walltimes <- df_walltimes |> tibble::as_tibble()
    # rm always col or powerset outcome without the "always" processes
    col_idx_keep <- grep(pattern = "always$", x = names(df_walltimes))
    df_walltimes <- df_walltimes |>
      dplyr::select(col_idx_keep)
    # bind n_inputs
    df_walltimes <- df_walltimes |>
      dplyr::mutate(n_inputs = 1:20)
    df_walltimes <- tidyr::gather(df_walltimes,
                                  key = "advanced_opts",
                                  value = "est_walltime",
                                  n_inputs)
    # sec to hrs
    df_walltimes <- df_walltimes |>
      dplyr::mutate(est_walltime = est_walltime / 3600)
    p <- ggplot2::ggplot(df_walltimes, ggplot2::aes(x = n_inputs, 
                                                    y = est_walltime, 
                                                    color = advanced_opts)) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "MolEvolvR estimated runtimes",
        x = "Number of inputs",
        y = "Estimated walltime (hours)"
      )
    return(p)
}
