############################################
##### Following code is from graph_ts: #####
############################################
args <- commandArgs(TRUE)
repo_path <- args[1]

# Get Conda env prefix from environment variable
conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")

# Set it as .libPaths so R uses packages from Conda env first
.libPaths(c(conda_lib, .libPaths()))
if (!requireNamespace("renv", quietly = TRUE)) {
  warning('R package renv not installed or not detected. If trimviz fails to find all R dependencies, please install renv, or run setup.sh and activate the trimViz2025_renv conda environment before calling trimviz.')
} else {
  renv::restore(project = repo_path) # should actually point to trimviz repo (not user's project dir)
}

# So there are 3 chances for loading dependencies:
# 1) renv (preferred)
# 2) Conda env (because conda_lib is prepended to .libPaths)
# 3) The user's default R libs

######################################
##### Now test for dependencies: #####
######################################

DEPENDENCIES <- c('ggplot2', 'ape', 'reshape2', 'gridExtra')
not_installed <- c()
for (dep in DEPENDENCIES) {
  if (!requireNamespace(dep, quietly = TRUE)) {
    not_installed <- c(not_installed, dep)
  }
}
if(length(not_installed) > 0) {
  msg <- paste0(not_installed, collapse=',')
  stop(paste0('Cannot find R dependencies: ', msg))
}
