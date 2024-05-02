# for now, we're using an env var, COMMON_SRC_ROOT, to specify this folder since
# the working directory is changed in many parts of the current molevolvr
# pipeline.
# to use this, construct paths like so: file.path(common_root, "path", "to", "file.R")
# for example, the reference for this file would be:
# file.path(common_root, "molevol_scripts", "R", "assign_job_queue.R")
.onLoad <- function(){
  common_root <- Sys.getenv("COMMON_SRC_ROOT")
}
