# ============================================================
# STANDALONE MICE PARALLELISM TEST  (run from project root)
# Loads data/analytic.rds, runs MICE (m=10), scores PREVENT on
# each imputation IN PARALLEL, and instruments worker overlap so
# we can SEE whether tasks ran concurrently or serially.
# Does NOT touch the rest of the pipeline or overwrite outputs.
# ============================================================
suppressPackageStartupMessages({
  library(survey); library(dplyr); library(mice)
  library(future); library(future.apply); library(preventr)
})
PROJ <- getwd()
t_start <- Sys.time()
cat("Loading analytic.rds...\n")
analytic <- readRDS(file.path(PROJ, "data", "analytic.rds"))

# ---- self-contained PREVENT wrapper (base model) -------------------------
prevent_risk <- function(age, sex, sbp, on_bpmed, total_chol, hdl_chol,
                         on_statin, diabetes, smoker, egfr, bmi = NULL,
                         model = "base") {
  winsorise <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  n <- length(age)
  bmi <- if (is.null(bmi)) rep(NA_real_, n) else bmi
  d <- data.frame(
    age=winsorise(age,30,79), sex=ifelse(sex=="Female","female","male"),
    sbp=winsorise(sbp,90,180), bp_tx=as.integer(ifelse(is.na(on_bpmed),FALSE,on_bpmed)),
    total_c=winsorise(total_chol,130,320), hdl_c=winsorise(hdl_chol,20,100),
    statin=as.integer(ifelse(is.na(on_statin),FALSE,on_statin)),
    dm=as.integer(ifelse(is.na(diabetes),FALSE,diabetes)),
    smoking=as.integer(ifelse(is.na(smoker),FALSE,smoker)),
    egfr=winsorise(egfr,15,150), bmi=winsorise(bmi,18.5,40), stringsAsFactors=FALSE)
  req <- c("age","sbp","total_c","hdl_c","egfr","bmi")
  ok <- stats::complete.cases(d[,req]) & age >= 30
  out <- rep(NA_real_, n); if (!any(ok)) return(out)
  dd <- d[ok,,drop=FALSE]; dd$preventr_id <- seq_len(nrow(dd))
  res <- preventr::estimate_risk(use_dat=dd, model=model, time="10yr",
                                 quiet=TRUE, collapse=TRUE, add_to_dat=FALSE, progress=FALSE)
  pos <- which(ok); out[pos[res$preventr_id]] <- as.numeric(res$total_cvd)*100
  out
}

# ---- MICE (uses futuremice internally; that's fine) -----------------------
M <- 10
mice_vars <- c("total_chol","hdl_chol","hba1c","egfr","bmi")
predictors <- c("age","sex","sbp_mean","diabetes_dx","smoking_status",
                "uacr","race_eth","education","pir")
mice_in <- analytic %>%
  select(all_of(unique(c(mice_vars, predictors,"sdmvpsu","sdmvstra",
                         "mec_wt_pooled","bp_nonattained","on_htn_meds","ckm_stage_f")))) %>%
  mutate(across(where(is.character), as.factor))

cat("Running MICE (m=10)...\n")
t_mice <- system.time(
  imp <- mice::futuremice(mice_in, m=M, method="pmm",
                          n.core=min(10, parallel::detectCores()-1),
                          parallelseed=2025)
)
cat(sprintf("MICE done in %.0fs\n", t_mice["elapsed"]))

# IMPORTANT: reset the plan AFTER futuremice, THEN set our own for scoring.
future::plan(future::sequential)               # clear futuremice's plan
# Trim each completed dataset to ONLY scoring columns (the key fix: avoids
# exporting full 17k-row x 17-col frames x10 to every worker).
score_cols <- c("age","sex","sbp_mean","on_htn_meds","total_chol","hdl_chol",
                "diabetes_dx","smoking_status","egfr","bmi",
                "sdmvpsu","sdmvstra","mec_wt_pooled")
completed_list <- lapply(seq_len(M), function(k) {
  ci <- mice::complete(imp, k)
  ci[, intersect(score_cols, names(ci)), drop = FALSE]
})

PROG <- file.path(PROJ, "logs", "test_mice_progress.txt")
cat(sprintf("TEST start %s\n", format(Sys.time())), file=PROG, append=FALSE)

score_one <- function(k, ci, src_fn, prog) {
  t0 <- Sys.time()
  ps <- src_fn(ci$age, ci$sex, ci$sbp_mean, ci$on_htn_meds %in% TRUE,
               ci$total_chol, ci$hdl_chol, FALSE,
               ci$diabetes_dx %in% TRUE, ci$smoking_status=="Current",
               ci$egfr, bmi=ci$bmi)
  di <- survey::svydesign(ids=~sdmvpsu, strata=~sdmvstra, weights=~mec_wt_pooled,
                          data=cbind(ci, ps=ps), nest=TRUE)
  m <- survey::svymean(~I(ps>=20), subset(di, !is.na(ps)), na.rm=TRUE)
  cat(sprintf("[%s] worker pid=%d imp=%d started~%s done\n",
              format(Sys.time(),"%H:%M:%OS2"), Sys.getpid(), k,
              format(t0,"%H:%M:%OS2")), file=prog, append=TRUE)
  list(est=as.numeric(m)[1], se=as.numeric(survey::SE(m))[1],
       pid=Sys.getpid(), start=as.numeric(t0), end=as.numeric(Sys.time()))
}

cat("Setting multisession plan for scoring...\n")
future::plan(future::multisession, workers=min(10, parallel::detectCores()-1))
cat("nbrOfWorkers =", nbrOfWorkers(), "\n")

t_score <- system.time(
  results <- future.apply::future_Map(
    function(k, ci) score_one(k, ci, prevent_risk, PROG),
    seq_len(M), completed_list,
    future.seed=TRUE, future.packages=c("survey","preventr"))
)
future::plan(future::sequential)

# ---- overlap diagnosis ----------------------------------------------------
pids   <- sapply(results, `[[`, "pid")
starts <- sapply(results, `[[`, "start")
ends   <- sapply(results, `[[`, "end")
span   <- max(ends) - min(starts)         # wall-clock for all tasks
sumdur <- sum(ends - starts)              # summed per-task durations
cat("\n================ PARALLELISM DIAGNOSIS ================\n")
cat(sprintf("distinct worker PIDs : %d (of %d tasks)\n", length(unique(pids)), M))
cat(sprintf("scoring wall-clock   : %.0fs\n", t_score["elapsed"]))
cat(sprintf("sum of task durations: %.0fs\n", sumdur))
cat(sprintf("overlap ratio        : %.1fx  (>1 => parallel; ~1 => serial)\n",
            sumdur / span))
cat(sprintf("VERDICT: %s\n", ifelse(sumdur/span > 1.5, "PARALLEL ✓", "SERIAL ✗")))
cat(sprintf("total test time      : %.0fs\n",
            as.numeric(difftime(Sys.time(), t_start, units="secs"))))
cat("======================================================\n")