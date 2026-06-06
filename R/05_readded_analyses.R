# =============================================================================
# SEGMENT 5 — RE-ADDED ANALYSES
# Brings back the analyses that exist in the reference pipeline but were dropped:
#   5A. Race-stratified primary regression (Model 3 within each race group)
#   5B. Race x PIR interaction test     -- Bonferroni family of 3,
#   5C. Race x PREVENT interaction test -- alpha = 0.05/3 = 0.0167
#       (UACR x CKM interaction from Segment 4 is the third member; its
#        significance verdict is re-evaluated here at the shared alpha.)
#   5D. UACR x CKM cell-level attainment + stratified PRs
#   5E. Antihypertensive drug-class distribution (G-J only; L excluded)
# All survey-weighted; modified Poisson for PRs; NCHS suppression applied.
# =============================================================================

stopifnot(exists("log_step"))
suppressPackageStartupMessages({
  library(survey); library(dplyr); library(tidyr); library(purrr); library(readr)
})
options(survey.lonely.psu = "adjust")

analytic <- get("analytic", envir = globalenv())
OUT <- file.path(PROJ_ROOT, "output")
BONF_ALPHA <- 0.05 / 3   # 3 pre-specified interaction tests

# Rebuild design objects local to this segment (kept independent of Segment 4)
design <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~mec_wt_pooled,
                    data = analytic, nest = TRUE)
design_m3 <- subset(design, nhanes_cycle != "L" & !is.na(antihtn_drug_classes))

fmt_p <- function(p) formatC(p, format = "e", digits = 10)
nchs_flag <- function(estimate, se, n_unwt) {
  rse <- ifelse(estimate != 0, 100 * se / abs(estimate), NA_real_)
  (!is.na(n_unwt) & n_unwt < 30) | (!is.na(rse) & rse > 30)
}
pr_table <- function(model) {
  s <- summary(model); co <- s$coefficients
  est <- co[, "Estimate"]; se <- co[, "Std. Error"]; p <- co[, ncol(co)]
  tibble(term = rownames(co), PR = exp(est),
         ci_lo = exp(est - 1.96*se), ci_hi = exp(est + 1.96*se),
         SE_logPR = se, p = p, p_fmt = fmt_p(p))
}
write_tbl <- function(tbl, name) {
  readr::write_csv(tbl, file.path(OUT, paste0(name, ".csv")))
  log_ok(sprintf("  wrote %s.csv (%d rows)", name, nrow(tbl)))
}

# Full Model-3 RHS including the new objective measures
rhs_full <- paste(
  "age10 + sex_f + race_eth_f + pir + education_f + uninsured +",
  "bmi5 + smoking_f + alcohol_f + egfr_cat_f + uacr_cat_f + hba1c +",
  "antihtn_drug_classes + ckm_stage_f + diet_sodium_500 +",
  "cotinine_log + pa_met_min_wk_1000"
)

# =====================================================================
# 5A — Race-stratified regression (Model 3 within each race group)
# Drops race_eth_f from RHS (constant within stratum).
# =====================================================================
log_step("Analysis 5A: Race-stratified Model 3")
rhs_norace <- gsub("\\+ race_eth_f ", "", rhs_full)
race_levels <- levels(analytic$race_eth_f)
if (exists("progress_init")) progress_init(length(race_levels), "Race strata")
strat_rows <- map_dfr(race_levels, function(rl) {
  on.exit(if (exists("progress_tick")) progress_tick(sprintf("(%s)", rl)), add = TRUE)
  d_r <- subset(design_m3, race_eth_f == rl)
  n_r <- nrow(d_r$variables)
  if (n_r < 30) {
    log_warn(sprintf("  race='%s' n=%d (<30) — suppressed", rl, n_r))
    return(tibble(term = NA, PR = NA, ci_lo = NA, ci_hi = NA, SE_logPR = NA,
                  p = NA, p_fmt = NA, race_eth = rl, n = n_r, suppressed = TRUE))
  }
  m <- tryCatch(
    svyglm(as.formula(paste("bp_nonattained ~", rhs_norace)),
           design = d_r, family = quasipoisson(link = "log")),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  pr_table(m) %>% mutate(race_eth = rl, n = n_r, suppressed = FALSE)
})
write_tbl(strat_rows, "analysis5A_race_stratified_regression")

# =====================================================================
# 5B — Race x PIR interaction (Bonferroni member)
# =====================================================================
log_step("Analysis 5B: Race x PIR interaction")
f_racepir <- as.formula(paste("bp_nonattained ~", rhs_full, "+ race_eth_f:pir"))
m_rp <- svyglm(f_racepir, design = design_m3, family = quasipoisson(link = "log"))
p_rp <- regTermTest(m_rp, ~ race_eth_f:pir)$p
write_tbl(tibble(
  test = "race_eth x PIR", p_interaction = as.numeric(p_rp),
  p_interaction_full = fmt_p(as.numeric(p_rp)),
  bonferroni_alpha = BONF_ALPHA, significant = as.numeric(p_rp) < BONF_ALPHA
), "analysis5B_race_pir_interaction")

# =====================================================================
# 5C — Race x PREVENT interaction (Bonferroni member)
# PREVENT score is continuous; restrict to rows with non-missing prevent_score.
# =====================================================================
log_step("Analysis 5C: Race x PREVENT interaction")
design_prev <- subset(design_m3, !is.na(prevent_score))
design_prev <- update(design_prev, prevent10 = prevent_score / 10)
f_racepv <- as.formula(paste("bp_nonattained ~", rhs_full, "+ race_eth_f:prevent10"))
m_rv <- tryCatch(
  svyglm(f_racepv, design = design_prev, family = quasipoisson(link = "log")),
  error = function(e) { log_warn(paste("  5C model issue:", conditionMessage(e))); NULL })
if (!is.null(m_rv)) {
  p_rv <- regTermTest(m_rv, ~ race_eth_f:prevent10)$p
  write_tbl(tibble(
    test = "race_eth x PREVENT", p_interaction = as.numeric(p_rv),
    p_interaction_full = fmt_p(as.numeric(p_rv)),
    bonferroni_alpha = BONF_ALPHA, significant = as.numeric(p_rv) < BONF_ALPHA
  ), "analysis5C_race_prevent_interaction")
}

# Re-state the UACR x CKM verdict at the shared Bonferroni alpha (third member)
uacr_ckm_path <- file.path(OUT, "analysis5_uacr_ckm_interaction.csv")
if (file.exists(uacr_ckm_path)) {
  uc <- readr::read_csv(uacr_ckm_path, show_col_types = FALSE)
  uc$bonferroni_alpha <- BONF_ALPHA
  uc$significant <- uc$p_interaction < BONF_ALPHA
  readr::write_csv(uc, uacr_ckm_path)
  log_ok("  re-evaluated UACR x CKM verdict at shared Bonferroni alpha (0.0167)")
}

# =====================================================================
# 5D — UACR x CKM cell-level attainment + stratified PRs
# =====================================================================
log_step("Analysis 5D: UACR x CKM cells & stratified PRs")
cell_means <- svyby(~bp_attained, ~ckm_stage_f + uacr_cat_f, design,
                    svymean, na.rm = TRUE)
cell_n <- analytic %>% count(ckm_stage_f, uacr_cat_f, name = "n_unwt")
cell_tbl <- as_tibble(cell_means) %>%
  rename(pct_attained = bp_attainedTRUE, se = se.bp_attainedTRUE) %>%
  left_join(cell_n, by = c("ckm_stage_f","uacr_cat_f")) %>%
  mutate(pct_nonatt = 1 - pct_attained,
         ci_lo = pct_attained - 1.96*se, ci_hi = pct_attained + 1.96*se,
         suppressed = nchs_flag(pct_attained, se, n_unwt))
write_tbl(cell_tbl, "analysis5D_uacr_ckm_prevalence")

# Stratified PR: within each CKM stage, PR of non-attainment for A2/A3 vs A1
strat_pr <- map_dfr(levels(analytic$ckm_stage_f), function(st) {
  d_s <- subset(design, ckm_stage_f == st)
  n_s <- nrow(d_s$variables)
  m <- tryCatch(svyglm(bp_nonattained ~ uacr_cat_f + age10 + sex_f,
                       design = d_s, family = quasipoisson(link="log")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  pt <- pr_table(m)
  pt <- pt[grepl("uacr_cat_f", pt$term), ]
  pt %>% transmute(ckm_stage = st, n_stage = n_s,
                   uacr_vs_a1 = gsub("uacr_cat_f", "", term),
                   PR, ci_lo, ci_hi, p, p_fmt)
})
write_tbl(strat_pr, "analysis5D_uacr_ckm_stratified_pr")

# =====================================================================
# 5E — Antihypertensive drug-class distribution (G-J only; L excluded)
# =====================================================================
log_step("Analysis 5E: Drug-class distribution (G-J)")
d_gj <- subset(design, nhanes_cycle != "L" & !is.na(antihtn_drug_classes))
# Distribution among treated participants
d_treated <- subset(d_gj, on_htn_meds %in% TRUE)
class_levels <- update(d_treated,
                       class_cat = factor(pmin(antihtn_drug_classes, 4),
                                          levels = 0:4, labels = c("0","1","2","3","4+")))
dist <- svymean(~class_cat, class_levels, na.rm = TRUE)
dist_n <- d_treated$variables %>%
  mutate(class_cat = pmin(antihtn_drug_classes, 4)) %>%
  count(class_cat, name = "n_unwt")
dist_tbl <- tibble(
  class_count = names(dist),
  pct = as.numeric(dist) * 100,
  se = as.numeric(SE(dist)) * 100
) %>%
  mutate(ci_lo = pct - 1.96*se, ci_hi = pct + 1.96*se)
write_tbl(dist_tbl, "analysis5E_drug_class_distribution")

# aTRH prevalence among treated (uses SAP definition already in analytic$atrh)
atrh_design <- subset(d_gj, on_htn_meds %in% TRUE & !is.na(atrh))
atrh_est <- svymean(~atrh, atrh_design, na.rm = TRUE)
write_tbl(tibble(
  metric = "aTRH prevalence among treated (G-J)",
  prop = as.numeric(atrh_est)[1], se = as.numeric(SE(atrh_est))[1]
), "analysis5E_atrh_prevalence")

# =====================================================================
# 5F — PREVENT pharmacotherapy-eligibility thresholds (7.5% / 10%)
# DISTINCT from the 20% Stage-3 staging cutoff: these are TREATMENT-decision
# thresholds (>=7.5% statin consideration; >=10% antihypertensive/risk-factor
# intensification per 2025 ACC/AHA). Denominator matches the reference pipeline:
#   Stage 2 + BP non-attained + age 30-79 + no clinical CVD (PREVENT-calculable,
#   treatment-relevant population; figure N ~ 5,455).
# =====================================================================
log_step("Analysis 5F: PREVENT pharmacotherapy thresholds (7.5/10%)")
rx_design <- subset(
  design,
  ckm_stage == 2 &
    bp_nonattained == 1 &
    age >= 30 & age <= 79 &
    !(cvd_history %in% TRUE) &
    !is.na(prevent_score)
)
rx_n <- nrow(rx_design$variables)
log_step(sprintf("  pharmacotherapy denominator N = %d", rx_n))

rx_tbl <- purrr::map_dfr(c(7.5, 10), function(th) {
  dd <- update(rx_design, exceed = as.integer(prevent_score >= th))
  m  <- svymean(~exceed, dd, na.rm = TRUE)
  est <- as.numeric(m)[1]; se <- as.numeric(SE(m))[1]
  tibble(
    threshold = th,
    pct_above = 100 * est,
    se = 100 * se,
    ci_lo = 100 * (est - 1.96 * se),
    ci_hi = 100 * (est + 1.96 * se),
    n_denominator = rx_n,
    label = sprintf("%.1f%% qualify for pharmacotherapy intensification", 100 * est),
    suppressed = nchs_flag(est, se, rx_n)
  )
})
write_tbl(rx_tbl, "analysis5F_pharmacotherapy_thresholds")

# Also stratified by race (mirrors reference analysis6_prevent_by_race)
rx_by_race <- purrr::map_dfr(levels(analytic$race_eth_f), function(rl) {
  d_r <- subset(rx_design, race_eth_f == rl)
  n_r <- nrow(d_r$variables)
  if (n_r < 30) return(tibble(race_eth = rl, n = n_r, threshold = NA,
                              pct_above = NA, ci_lo = NA, ci_hi = NA,
                              suppressed = TRUE))
  purrr::map_dfr(c(7.5, 10), function(th) {
    dd <- update(d_r, exceed = as.integer(prevent_score >= th))
    m  <- svymean(~exceed, dd, na.rm = TRUE)
    est <- as.numeric(m)[1]; se <- as.numeric(SE(m))[1]
    tibble(race_eth = rl, n = n_r, threshold = th,
           pct_above = 100 * est, ci_lo = 100*(est-1.96*se),
           ci_hi = 100*(est+1.96*se), suppressed = nchs_flag(est, se, n_r))
  })
})
write_tbl(rx_by_race, "analysis5F_pharmacotherapy_by_race")

log_ok("Segment 5 complete: re-added analyses written to /output.")