
# example: list_opts2procs <- make_opts2procs
make_opts2procs <- function() {
  opts2processes <- list(
    "homology_search" = c("dblast", "dblast_cleanup"),
    "domain_architecture" = c("iprscan", "ipr2lineage", "ipr2da"),
    "any" = c("blast_clust", "clust2table") # processes always present agnostic of input opts 
  )
  return(opts2processes)
}

# example:
# in_opts <- c("homology_search", "domain_architecture")
# procs <- map_input_opts2procs(in_opts)
map_input_opts2procs <- function(input_opts) {
  list_opts2proc <- make_opts2procs()
  idx <- which(input_opts %in% names(list_opts2proc))
  return(list_opts2proc[idx] |> unlist())
}

in_opts <- c("homology_search", "domain_architecture")
procs <- map_input_opts2procs(in_opts)

### example code to get the median process walltimes from logs
# run this function outside of container
get_proc_medians <- function() {
  source("/data/molevolvr_transfer/molevolvr_dev/molevol_scripts/R/metrics.R")
  # aggregate logs from
  path_prod_results <- "/data/molevolvr_transfer/hpc-cluster-tests/job_results"
  path_log_data <- "/data/molevolvr_transfer/molevolvr_dev/molevol_scripts/log_data/prod_logs.rda"
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
  ?format.POSIXct()
}

# example: write_proc_medians_table()
write_proc_medians_table <- function(
    filepath = "/data/molevolvr_transfer/molevolvr_dev/molevol_scripts/log_data/prod_process_medians.tsv"
) {

proc_medians <- get_proc_medians()
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

# molevolvr backend process weight list based on the median walltimes 
# for each process (proc)
# example: get_proc_weights()
get_proc_weights <- function() {
  proc_weights <- list(
    "dblast" = 2810,
    "iprscan" = 1016,
    "dblast_cleanup" = 79,
    "ipr2lineage" = 18,
    "ipr2da" = 12,
    "blast_clust" = 2,
    "clust2table" = 2
  )
  return(proc_weights)
}

# example: .test_input_opts_to_est_walltime("homology_search")
.test_input_opts_to_est_walltime <- function(input_opts, n_inputs = 1L) {
  proc_weights <- get_proc_weights()
  # sort the list by names and convert to vec
  proc_weights <- proc_weights[order(names(proc_weights))] |> unlist()
  all_procs <- names(proc_weights) |> sort()
  procs_from_opts <- map_input_opts2procs(input_opts)
  procs_from_opts <- sort(procs_from_opts)
  binary_proc_vec  <- dplyr::if_else(all_procs %in% procs_from_opts, 1L, 0L)
  # dot product of weights and procs to run multiplied by the number of inputs
  est_walltime <- n_inputs * (binary_proc_vec %*% proc_weights) |> as.numeric()
  return(est_walltime)
}
