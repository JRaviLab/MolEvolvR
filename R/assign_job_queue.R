
# for now, we're using an env var, COMMON_SRC_ROOT, to specify this folder since
# the working directory is changed in many parts of the current molevolvr
# pipeline.
# to use this, construct paths like so: file.path(common_root, "path", "to", "file.R")
# for example, the reference for this file would be:
# file.path(common_root, "molevol_scripts", "R", "assignJobQueue.R")
common_root <- Sys.getenv("COMMON_SRC_ROOT")

#' mapOption2Process
#' 
#' @description
#' Construct list where names (MolEvolvR advanced options) point to processes
#'
#' @importFrom rlang warn abort inform
#'
#' @return list where names (MolEvolvR advanced options) point to processes
#'
#' example: list_opts2procs <- mapOption2Process
#' @export
mapOption2Process <- function() {
  tryCatch({
    opts2processes <- list(
      "homology_search" = c("dblast", "dblast_cleanup"),
      "domain_architecture" = c("iprscan", "ipr2lineage", "ipr2da"),
      # processes always present agnostic of advanced options
      "always" = c("blast_clust", "clust2table")
    )
    return(opts2processes)
  }, error = function(e) {
    rlang::abort(paste("Error: ", e$message), class = "Opts_to_process_error")
  }, warning = function(w) {
    rlang::warn(paste("Warning: ", w$message),
                class = "Opts_to_process_warning")
  })

}

#' Use MolEvolvR advanced options to get associated processes
#'
#' @param advanced_opts character vector of MolEvolvR advanced options
#'
#' @importFrom rlang warn abort inform
#'
#' @return character vector of process names that will execute given
#' the advanced options
#'
#' example:
#' advanced_opts <- c("homology_search", "domain_architecture")
#' procs <- mapAdvOption2Process(advanced_opts)
#' @export
mapAdvOption2Process <- function(advanced_opts) {
  if (!is.character(advanced_opts)) {
    rlang::abort("Argument must be a character vector!",
                 class = "validation_error")
  }
  tryCatch({
    # append 'always' to add procs that always run
    advanced_opts <- c(advanced_opts, "always")
    opts2proc <- make_opts2procs()
    # setup index for opts2proc based on advanced options
    idx <- which(names(opts2proc) %in% advanced_opts)
    # extract processes that will run
    procs <- opts2proc[idx] |> unlist()
    return(procs)
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error: ", e$message),
      class = "map_advanced_opts2procs_error",
      advanced_opts = advanced_opts
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "map_advanced_opts2procs_warning",
      advanced_opts = advanced_opts
    )
  })

}

#' Scrape MolEvolvR logs and calculate median processes
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' directory
#'
#' @importFrom dplyr across everything select summarise
#' @importFrom rlang warn abort inform
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
  tryCatch({
    # Check if dir_job_results is a character string
    if (!is.character(dir_job_results) || length(dir_job_results) != 1) {
      rlang::abort("Input 'dir_job_results' must be a single character string.",
                   class = "validation_error")
    }

    # Check if dir_job_results exists
    if (!dir.exists(dir_job_results)) {
      rlang::abort(paste("The directory", dir_job_results, "does not exist."),
                         class = "file_error")
    }

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
  }, error = function(e) {
    rlang::abort(paste("Encountered an error: ", e$message),
                 class = "processing_error")
  }, warning = function(w) {
    rlang::warn(paste("Warning: ", w$message), class = "processing_warning")
  })

}

#' writeProcessRuntime2TSV
#' 
#' @description
#' Write a table of 2 columns: 1) process and 2) median seconds
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' @param filepath path to save tsv file
#'
#' @importFrom dplyr arrange desc everything
#' @importFrom tibble as_tibble
#' @importFrom readr write_tsv
#' @importFrom tidyr pivot_longer
#' @importFrom rlang warn abort inform
#'
#' @return [tbl_df] 2 columns: 1) process and 2) median seconds
#'
#' example: writeProcessRuntime2TSV(
#'   "/data/scratch/janani/molevolvr_out/",
#'   "/data/scratch/janani/molevolvr_out/log_tbl.tsv"
#' )
#' @export
writeProcessRuntime2TSV <- function(dir_job_results, filepath) {
  tryCatch({
    # Error handling for input arguments
    if (!is.character(dir_job_results) || length(dir_job_results) != 1) {
      rlang::abort("Input 'dir_job_results' must be a single character string.",
                   class = "validation_error")
    }

    if (!dir.exists(dir_job_results)) {
      rlang::abort(paste("The directory", dir_job_results, "does not exist."),
                   class = "file_error")
    }

    if (!is.character(filepath) || length(filepath) != 1) {
      rlang::abort("Input 'filepath' must be a single character string.",
                   class = "validation_error")
    }
    df_proc_medians <- get_proc_medians(dir_job_results) |>
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
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error: ", e$message),
      class = "processing_error",
      dir_job_results = dir_job_results,
      filepath = filepath
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "processing_warning",
      dir_job_results = dir_job_results,
      filepath = filepath
    )
  })

}

#' Compute median process runtimes, then write a YAML list of the processes and
#' their median runtimes in seconds to the path specified by 'filepath'.
#'
#' The default value of filepath is the value of the env var
#' MOLEVOLVR_PROC_WEIGHTS, which getProcessRuntimeWeights() also uses as its default
#' read location.
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results directory
#' @param filepath [chr] path to save YAML file; if NULL,
#'                 uses ./molevol_scripts/log_data/job_proc_weights.yml
#'
#' @importFrom yaml write_yaml
#' @importFrom rlang warn abort inform
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
  if (is.null(filepath)) {
    filepath <- file.path(common_root, "molevol_scripts", "log_data", "job_proc_weights.yml")
  }
  tryCatch({
    # Error handling for dir_job_results arguments
    if (!is.character(dir_job_results) || length(dir_job_results) != 1) {
      rlang::abort(
        message = "Input 'dir_job_results' must be a single character string.",
        class = "validation_error",
        dir_job_results = dir_job_results
      )
    }

    if (!dir.exists(dir_job_results)) {
      rlang::abort(
        message = paste("The directory", dir_job_results, "does not exist."),
        class = "file_error",
        dir_job_results = dir_job_results
      )
    }

    if (is.null(filepath)) {
      filepath <- file.path(common_root,
                            "molevol_scripts",
                            "log_data",
                            "job_proc_weights.yml")
    }
    if (!is.character(filepath) || length(filepath) != 1) {
      rlang::abort(
        message = "Input 'filepath' must be a single character string.",
        class = "validation_error",
        filepath = filepath
      )
    }

    medians <- calculateProcessRuntime(dir_job_results)
    yaml::write_yaml(medians, filepath)
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error: ", e$message),
      class = "processing_error",
      dir_job_results = dir_job_results,
      filepath = filepath
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "processing_warning",
      dir_job_results = dir_job_results,
      filepath = filepath
    )
  })
}

#' Quickly get the runtime weights for MolEvolvR backend processes
#'
#' @param dir_job_results [chr] path to MolEvolvR job_results
#' directory
#'
#' @importFrom stringr str_glue str_trim
#' @importFrom yaml read_yaml
#' @importFrom rlang warn abort inform
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
      rlang::abort(
        message = stringr::str_glue("medians_yml_path is empty
                                    ({medians_yml_path}), returning default weights"),
        class = "input_error",
        medians_yml_path = medians_yml_path
      )
    }

    proc_weights <- yaml::read_yaml(medians_yml_path)

    if (!is.list(proc_weights) || length(proc_weights) == 0) {
      rlang::abort(
        message = "The loaded YAML file does not
                  contain valid process weights.",
        class = "file_error",
        medians_yml_path = medians_yml_path
      )
    }
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

#' calculateEstimatedWallTimeFromOpts
#' 
#' @description
#' Given MolEvolvR advanced options and number of inputs,
#' calculate the total estimated walltime for the job
#'
#' @param advanced_opts character vector of MolEvolvR advanced options
#' (see mapOption2Process for the options)
#' @param n_inputs total number of input proteins
#'
#' @importFrom dplyr if_else
#' @importFrom stringr str_glue
#' @importFrom rlang warn abort inform
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

  tryCatch({
    # to calculate est walltime for a homology search job, the number of hits
    # must be provided
    validation_fail <- is.null(n_hits) && "homology_search" %in% advanced_opts
    stopifnot(!validation_fail)

    # Validate advanced_opts
    if (!is.character(advanced_opts)) {
      rlang::abort(
        message = "Argument 'advanced_opts' must be a character vector.",
        class = "validation_error",
        advanced_opts = advanced_opts
      )
    }

    # Validate n_inputs
    if (!is.numeric(n_inputs) || length(n_inputs) != 1 || n_inputs <= 0) {
      rlang::abort(
        message = "Argument 'n_inputs'
                  must be a single positive numeric value.",
        class = "validation_error",
        n_inputs = n_inputs
      )
    }

    # Validate n_hits if homology_search is in advanced_opts
    if ("homology_search" %in% advanced_opts &&
          (is.null(n_hits) || !is.numeric(n_hits) ||
             length(n_hits) != 1 || n_hits < 0)) {
      rlang::abort(
        message = "Argument 'n_hits' must be a single non-negative numeric
        value when 'homology_search' is in 'advanced_opts'.",
        class = "validation_error",
        n_hits = n_hits
      )
    }

  # Get process weights
    proc_weights <- writeProcessRuntime2YML()
    if (!is.list(proc_weights)) {
      rlang::abort(
        message = "Process weights could not be retrieved correctly.",
        class = "processing_error"
      )
    }

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
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error: ", e$message),
      class = "processing_error",
      advanced_opts = advanced_opts,
      n_inputs = n_inputs,
      n_hits = n_hits
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "processing_warning",
      advanced_opts = advanced_opts,
      n_inputs = n_inputs,
      n_hits = n_hits
    )
  })

}


#' Decision function to assign job queue
#'
#' @param t_sec_estimate estimated number of seconds a job will process
#' (from calculateEstimatedWallTimeFromOpts	())
#' @param t_long threshold value that defines the lower bound for assigning a
#' job to the "long queue"
#'
#' @importFrom rlang warn abort inform
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
  tryCatch({
    # Validate t_sec_estimate
    if (!is.numeric(t_sec_estimate) || length(t_sec_estimate) != 1) {
      rlang::abort(
        message = "Argument 't_sec_estimate' must be a single numeric value.",
        class = "validation_error",
        t_sec_estimate = t_sec_estimate
      )
    }

    # Validate t_cutoff
    if (!is.numeric(t_cutoff) || length(t_cutoff) != 1 || t_cutoff < 0) {
      rlang::abort(
        message = "Argument 't_cutoff' must be a
                  single non-negative numeric value.",
        class = "validation_error",
        t_cutoff = t_cutoff
      )
    }


    queue <- ifelse(t_sec_estimate > t_cutoff, "long", "short")
    return(queue)
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error: ", e$message),
      class = "processing_error",
      t_sec_estimate = t_sec_estimate,
      t_cutoff = t_cutoff
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "processing_warning",
      t_sec_estimate = t_sec_estimate,
      t_cutoff = t_cutoff
    )
  })

}

#' Plot the estimated runtimes for different advanced options and number
#' of inputs
#'
#' this function was just for fun; very, very messy code
#'
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 aes geom_line ggplot labs
#' @importFrom tibble as_tibble
#' @importFrom rlang warn abort inform
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
  }, error = function(e) {
    rlang::abort(
      message = paste("Encountered an error:", e$message),
      .internal = TRUE
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning:", w$message),
      .internal = TRUE
    )
  })

}
