# =============================================================================
# SEGMENT 1 — ENVIRONMENT & PACKAGE SETUP
# Initialises renv in the project folder, installs all required packages into
# the project-local library, and snapshots versions (incl. R version) to
# renv.lock. Records the RStudio version if discoverable.
# =============================================================================

stopifnot(exists("log_step"))  # must be sourced by run_study.R

# ---- Required packages ------------------------------------------------------
REQUIRED_PKGS <- c(
  "renv",        # environment capture
  "nhanesA",     # NHANES retrieval
  "preventr",    # official AHA PREVENT equations (Khan 2023)
  "survey",      # complex survey design + svyglm / svyvar / svyby
  "dplyr",       # data wrangling
  "tidyr",       # reshaping
  "purrr",       # functional iteration
  "stringr",     # string matching (drug names)
  "forcats",     # factor handling
  "mice",        # multiple imputation
  "car",         # VIF
  "parallel",    # core detection
  "future",      # parallel backend
  "future.apply",# parallel apply
  "digest",      # hashing for PREVENT cache invalidation
  "tibble",
  "readr"        # writing output tables
)

log_step("Bootstrapping renv in project library")

# Install renv first (into default lib) if absent, then initialise project env.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialise renv if not already present. bare=TRUE so we control installs.
if (!file.exists(file.path(PROJ_ROOT, "renv.lock")) &&
    !dir.exists(file.path(PROJ_ROOT, "renv"))) {
  renv::init(bare = TRUE, restart = FALSE)
  log_ok("renv initialised (bare).")
} else {
  log_step("renv already present; activating.")
  source(file.path(PROJ_ROOT, "renv", "activate.R"))
}

# ---- Install any missing required packages into the project library ---------
installed <- rownames(installed.packages())
to_install <- setdiff(REQUIRED_PKGS, installed)
if (length(to_install)) {
  log_step(sprintf("Installing %d package(s): %s",
                   length(to_install), paste(to_install, collapse = ", ")))
  renv::install(to_install)
} else {
  log_ok("All required packages already installed.")
}

# ---- Load to confirm availability ------------------------------------------
invisible(lapply(setdiff(REQUIRED_PKGS, c("renv", "parallel")), function(p) {
  suppressPackageStartupMessages(require(p, character.only = TRUE)) ||
    stop(sprintf("Package failed to load: %s", p))
}))
log_ok("All packages loaded successfully.")

# ---- Detect cores for later parallel sections -------------------------------
N_CORES <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
log_step(sprintf("Detected usable cores for parallel work: %d", N_CORES))

# ---- Record RStudio version if running inside RStudio -----------------------
rstudio_ver <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) as.character(rstudioapi::getVersion()) else NA_character_
}, error = function(e) NA_character_)
writeLines(c(
  sprintf("R version: %s", R.version.string),
  sprintf("Platform: %s", R.version$platform),
  sprintf("RStudio version: %s", ifelse(is.na(rstudio_ver), "not detected / not in RStudio", rstudio_ver)),
  sprintf("renv snapshot time: %s", format(Sys.time()))
), file.path(PROJ_ROOT, "logs", "environment.txt"))

# ---- Snapshot the lockfile --------------------------------------------------
log_step("Writing renv.lock snapshot (package + R versions)")
renv::snapshot(prompt = FALSE)
log_ok("renv.lock written.")

# Expose constants used downstream.
assign("N_CORES", N_CORES, envir = globalenv())