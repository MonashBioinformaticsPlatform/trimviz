args = commandArgs(trailingOnly=TRUE)
cat('\n Running setup R script \n')
MODE      <- args[1]  # 'build' or 'test'
TV_PATH <- args[2]  # e.g. path/to/trimviz/ to avoid changing user's libraries
DEPENDENCIES <- c('ggplot2', 'ape', 'reshape2', 'gridExtra')
env_name <- Sys.getenv("CONDA_DEFAULT_ENV")
cat("\nActive conda environment name:", env_name, "\n")
if(env_name != 'trimViz2025_renv'){
 exit('Conda env is not trimViz2025_renv.')
}
conda_path <- Sys.getenv("CONDA_PREFIX")
.libPaths(  c(paste0(conda_path, '/lib/R/library/'), .libPaths())  )
cat('\n Libpaths: ', .libPaths(), '\n')

if (!requireNamespace("renv", quietly = TRUE)) {
  cat('\n renv should have been installed in the conda environment... have you activated it?\n')
  exit('renv not installed')
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(renv)

if(MODE == 'build'){  # build from the TV repo's renv.lock
  renv::restore(project = TV_PATH)
} else if(MODE == 'build_from_scratch') { # This is the code used to make the renv.lock file originally, which is included in the git repo. 
  setwd(TV_PATH)
  renv::init(project=TV_PATH, bare=T)
  for(dep in DEPENDENCIES){
    #if (!requireNamespace(dep, quietly = TRUE)) {
      #install.packages(dep, lib='./Rlib_renv')
    cat('\n installing: ', dep, '\n')
    renv::install(dep)
  }
  renv::snapshot() # project=TV_PATH) # setting project path again puts the lockfile in project_path/project_path!
  cat('\n Snapshotted renv \n')
} else if (MODE == 'test'){
  renv::restore(project = TV_PATH)
  library(ggplot2)
  library(ape)
  library(reshape2)
  library(gridExtra)
  cat('\n Restored renv.\n')
  print(sessionInfo())
} else {
 cat('\n Set MODE to either build or test.\n')
 cat('\n Arguments: 1) MODE 2) PROJECT_PATH 3) RLIB_PATH.  e.g. cd to trimviz repo dir, then:  Rscript setup_r_env.R build ./ ./Rlibs \n')
}
