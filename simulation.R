# If you want to rerun the whole simulation set this variable
# to `TRUE`
RERUN = TRUE

# If you just want to have a SHORT test simulation with just
# 10 repetitions, then set `TEST = TRUE`. The full simulation
# is done with 10000 repetitions when setting `TEST = FALSE`.
TEST = FALSE

# Install packages with respective version.
source(here::here("R/install-pkgs-versions.R"))

# But, even the full benchmark with 10000 repetitions is not
# too computational expensive and takes about 20 - 40 minutes
# on a Linux notebook with an
#
#        Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz
#
# processor and using 7 of the available 8 cores.
source(here::here("R/setup.R"))

# MAKE SURE THAT THE WD IS SET TO THE ROOT OF THE REPO!
setwd(here::here())

# We want to have relative paths here to ensure that the
# code also works on other machines. Therefore, it is
# important to set the working directory to the root of
# the repository.
BATCHTOOLS_DIR = "batchtools/"

## Batchtools
## ===========================================================

library(batchtools)

# Remove the batchtools directory to rerun the simulations.
if (dir.exists(BATCHTOOLS_DIR) && RERUN) unlink("batchtools", recursive = TRUE)

# If the batchtools directory exists, the simulations are
# continued:
if (dir.exists(BATCHTOOLS_DIR) && (! RERUN)) {

  loadRegistry(BATCHTOOLS_DIR, writeable = TRUE, work.dir = here::here())
  submitJobs(findNotDone())

# Otherwise, a new directory is created and the simulation
# is started:
} else {
  reg = makeExperimentRegistry(
    file.dir = BATCHTOOLS_DIR,
    packages = c("checkmate", "pROC"),
    source   = FILES,
    seed     = 31415)

  reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = parallel::detectCores())

  saveRegistry(reg)

  source("R/add-experiments.R")
  submitJobs()

  # Get the job status:
  getJobStatus()
}

