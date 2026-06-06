#!/usr/bin/env Rscript
# =============================================================================
# CKM SYNDROME × BLOOD PRESSURE GOAL ATTAINMENT — MASTER RUNNER
# Replication pipeline: 2025 ACC/AHA BP goal attainment across CKM stages 2-4
# Data: Pooled NHANES G,H,I,J (2011-2018) + L (2021-2023)
# -----------------------------------------------------------------------------
# This master script sources each segment in order. Each segment:
#   - announces its stage,
#   - writes a timestamped line to logs/run_log.txt on success or error,
#   - HALTS the whole pipeline on any error (stop-on-error contract).
# =============================================================================

# --- Resolve project root regardless of how the script is launched ----------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
PROJ_ROOT <- get_script_dir()
setwd(PROJ_ROOT)

dir.create("logs",   showWarnings = FALSE)
dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

LOG_FILE <- file.path(PROJ_ROOT, "logs", "run_log.txt")

# --- Logging helpers ---------------------------------------------------------
log_msg <- function(level, msg) {
  line <- sprintf("[%s] %-7s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  cat(line, "\n")
  flush.console()
}
log_step  <- function(msg) log_msg("STEP",  msg)
log_ok    <- function(msg) log_msg("OK",    msg)
log_warn  <- function(msg) log_msg("WARN",  msg)
log_error <- function(msg) log_msg("ERROR", msg)

# ---- Console progress helpers (console-only; not written to run_log) --------
# Detect whether we're in an interactive/RStudio console (where \r overwrite is
# unreliable) vs a real terminal via Rscript. In the former, print periodic new
# lines instead of overwriting, so output never gets mangled.
.PROG_INTERACTIVE <- interactive()
.PROG_STATE <- new.env(parent = emptyenv())

# progress_init(total, label): start a progress run of `total` ticks.
progress_init <- function(total, label) {
  .PROG_STATE$total <- total
  .PROG_STATE$label <- label
  .PROG_STATE$n     <- 0L
  .PROG_STATE$t0    <- Sys.time()
  cat(sprintf("    %s: 0/%d ...\n", label, total))
  flush.console()
}

# progress_tick(extra): advance one step; refresh the status line.
progress_tick <- function(extra = "") {
  s <- .PROG_STATE
  s$n <- s$n + 1L
  pct <- round(100 * s$n / s$total)
  el  <- as.numeric(difftime(Sys.time(), s$t0, units = "secs"))
  # crude ETA from average pace so far
  eta <- if (s$n > 0) el / s$n * (s$total - s$n) else NA
  nb  <- round(20 * s$n / s$total)
  bar <- paste0(strrep("#", nb), strrep("-", 20 - nb))
  line <- sprintf("    %s [%s] %d/%d (%d%%) | %.0fs elapsed, ~%.0fs left %s",
                  s$label, bar, s$n, s$total, pct, el,
                  ifelse(is.na(eta), 0, eta), extra)
  if (.PROG_INTERACTIVE) {
    # RStudio: print a fresh line only at each 10% boundary to avoid spam
    if (s$n == s$total || pct %% 10 == 0) cat(line, "\n")
  } else {
    cat("\r", line, sep = ""); if (s$n == s$total) cat("\n")
  }
  flush.console()
}

# notice(msg): a one-off heads-up for long un-loopable operations.
notice <- function(msg) { cat(sprintf("    >> %s\n", msg)); flush.console() }

banner <- function(txt) {
  bar <- paste(rep("=", 78), collapse = "")
  cat("\n", bar, "\n  ", txt, "\n", bar, "\n", sep = "")
}

# Wrap a segment so any error is logged and then re-thrown to halt the run.
.SEG_TOTAL <- 5L
.SEG_N <- 0L
run_segment <- function(seg_path, seg_label) {
  .SEG_N <<- .SEG_N + 1L
  banner(sprintf("[Segment %d/%d] %s", .SEG_N, .SEG_TOTAL, seg_label))
  log_step(sprintf("BEGIN [%d/%d] %s", .SEG_N, .SEG_TOTAL, seg_label))
  t0 <- Sys.time()
  res <- tryCatch(
    {
      sys.source(seg_path, envir = globalenv())
      TRUE
    },
    error = function(e) {
      log_error(sprintf("%s FAILED: %s", seg_label, conditionMessage(e)))
      e
    }
  )
  if (inherits(res, "error")) {
    log_error("Pipeline halted. See message above. Fix and re-run.")
    quit(save = "no", status = 1L)
  }
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_ok(sprintf("END   [%d/%d] %s (%.1fs)", .SEG_N, .SEG_TOTAL, seg_label, dt))
}

# --- Fresh log header --------------------------------------------------------
cat(sprintf("CKM STUDY RUN LOG — started %s\n", format(Sys.time())),
    file = LOG_FILE, append = FALSE)
log_step(sprintf("Project root: %s", PROJ_ROOT))

# --- Run order ---------------------------------------------------------------
run_segment(file.path("R", "01_setup_packages.R"),   "SEGMENT 1: Environment & package setup (renv)")
run_segment(file.path("R", "02_load_merge.R"),        "SEGMENT 2: Download & merge NHANES variables")
run_segment(file.path("R", "03_derive_stage.R"),      "SEGMENT 3: Derive variables & assign CKM stages")
run_segment(file.path("R", "04_analyses.R"),          "SEGMENT 4: Survey-weighted analyses (1-10)")
run_segment(file.path("R", "05_readded_analyses.R"),  "SEGMENT 5: Re-added analyses (race, interactions, drug classes)")

banner("PIPELINE COMPLETE")
log_ok("All segments completed successfully. Outputs in /output, log in /logs.")