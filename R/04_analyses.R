# =============================================================================
# SEGMENT 4 — SURVEY-WEIGHTED ANALYSES (1-10)
# Builds survey design objects and runs the full statistical analysis plan.
# All estimates are survey-weighted. Modified Poisson (quasipoisson log link,
# robust variance) is the primary regression engine -> Prevalence Ratios.
# Output tables carry full (untruncated) p-values, SEs, and SMDs, and apply
# NCHS suppression flags (RSE > 30% or unweighted n < 30).
# =============================================================================

stopifnot(exists("log_step"))
suppressPackageStartupMessages({
  library(survey); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(forcats); library(car); library(mice)
})
options(survey.lonely.psu = "adjust")

analytic <- get("analytic", envir = globalenv())
OUT <- file.path(PROJ_ROOT, "output")

# ---- Survey design objects --------------------------------------------------
log_step("Constructing survey design objects")
design <- svydesign(
  ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~mec_wt_pooled,
  data = analytic, nest = TRUE
)
# Regression design (drop rows missing core model fields handled per-model via subset)
design_reg <- design

# Fasting subsample design (WTSAF2YR/5)
fasting <- analytic %>% filter(is_fasting, !is.na(saf_wt_pooled), saf_wt_pooled > 0)
design_fast <- if (nrow(fasting) > 30) {
  svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~saf_wt_pooled,
            data = fasting, nest = TRUE)
} else NULL
log_ok(sprintf("Designs built. Main n=%d; fasting n=%d",
               nrow(analytic), nrow(fasting)))

# ---- Shared helpers ---------------------------------------------------------
# NCHS suppression: flag if unweighted n < 30 OR RSE > 30%.
nchs_flag <- function(estimate, se, n_unwt) {
  rse <- ifelse(estimate != 0, 100 * se / abs(estimate), NA_real_)
  (!is.na(n_unwt) & n_unwt < 30) | (!is.na(rse) & rse > 30)
}

fmt_p <- function(p) formatC(p, format = "e", digits = 10)  # full, untruncated

# Modified Poisson PR extractor with robust CI
pr_table <- function(model) {
  s <- summary(model)
  co <- s$coefficients
  est <- co[, "Estimate"]; se <- co[, "Std. Error"]
  p   <- co[, ncol(co)]
  tibble(
    term = rownames(co),
    PR   = exp(est),
    CI_low = exp(est - 1.96 * se),
    CI_high= exp(est + 1.96 * se),
    SE_logPR = se,
    p_value = p,
    p_value_full = fmt_p(p)
  )
}

write_tbl <- function(tbl, name) {
  path <- file.path(OUT, paste0(name, ".csv"))
  readr::write_csv(tbl, path)
  log_ok(sprintf("  wrote %s (%d rows)", basename(path), nrow(tbl)))
}

# =====================================================================
# ANALYSIS 1 — Descriptive statistics + primary outcome
# =====================================================================
log_step("Analysis 1: Descriptive statistics & primary outcome")

# Primary outcome: proportion BP-attained, overall and by CKM stage
prop_overall <- svymean(~bp_attained, design, na.rm = TRUE)
by_stage <- svyby(~bp_attained, ~ckm_stage_f, design, svymean, na.rm = TRUE)
by_stage_n <- analytic %>% count(ckm_stage_f, name = "n_unwt")

primary_tbl <- as_tibble(by_stage) %>%
  rename(prop_attained = bp_attainedTRUE, se = se.bp_attainedTRUE) %>%
  left_join(by_stage_n, by = "ckm_stage_f") %>%
  mutate(
    ci_low = prop_attained - 1.96 * se,
    ci_high = prop_attained + 1.96 * se,
    suppress = nchs_flag(prop_attained, se, n_unwt)
  )

# P-trend: survey logistic with stage as continuous (2,3,4)
analytic$ckm_stage_num <- as.numeric(as.character(analytic$ckm_stage_f))
design_trend_stage <- update(design, ckm_stage_num = as.numeric(as.character(ckm_stage_f)))
trend_mod <- svyglm(bp_attained ~ ckm_stage_num, design = design_trend_stage,
                    family = quasibinomial())
p_trend <- summary(trend_mod)$coefficients["ckm_stage_num", 4]
primary_tbl$p_trend_across_stages <- fmt_p(p_trend)
write_tbl(primary_tbl, "analysis1_primary_bp_attainment_by_stage")

# Secondary: phenotype distribution among Stage 2+ (treated phenotypes)
pheno_design <- subset(design, bp_phenotype %in%
                         c("Controlled-Treated","Uncontrolled-Treated","Uncontrolled-Untreated"))
pheno_overall <- svymean(~factor(bp_phenotype), pheno_design, na.rm = TRUE)
pheno_by_stage <- svyby(~factor(bp_phenotype), ~ckm_stage_f, pheno_design,
                        svymean, na.rm = TRUE)
write_tbl(as_tibble(pheno_by_stage), "analysis1b_phenotype_by_stage")

# Descriptive means for continuous covariates by stage
cont_vars <- c("age","bmi","waist_cm","egfr","hba1c","uacr","sbp_mean",
               "dbp_mean","total_chol","hdl_chol","prevent_score","diet_sodium")
desc_rows <- map_dfr(cont_vars, function(v) {
  f <- as.formula(paste0("~", v))
  m <- tryCatch(svyby(f, ~ckm_stage_f, design, svymean, na.rm = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  as_tibble(m) %>% mutate(variable = v) %>%
    rename(mean = 2, se = 3)
})
write_tbl(desc_rows, "analysis1c_continuous_means_by_stage")

# =====================================================================
# ANALYSIS 2 — Phenotype table with survey-weighted SMDs
# =====================================================================
log_step("Analysis 2: SMD table by BP phenotype")

pheno_levels <- c("Controlled-Treated","Uncontrolled-Treated","Uncontrolled-Untreated")
pd <- subset(design, bp_phenotype %in% pheno_levels)

# Continuous SMD using svyvar() population variance (NOT SE^2)
smd_continuous <- function(var) {
  f <- as.formula(paste0("~", var))
  means <- svyby(f, ~bp_phenotype, pd, svymean, na.rm = TRUE)
  vars  <- svyby(f, ~bp_phenotype, pd, svyvar,  na.rm = TRUE)  # population variance
  mm <- setNames(means[[2]], means$bp_phenotype)
  vv <- setNames(vars[[2]],  vars$bp_phenotype)
  combos <- combn(pheno_levels, 2)
  smds <- apply(combos, 2, function(pr) {
    a <- pr[1]; b <- pr[2]
    abs(mm[a] - mm[b]) / sqrt((vv[a] + vv[b]) / 2)
  })
  tibble(variable = var, type = "continuous",
         max_smd = max(smds, na.rm = TRUE),
         detail = paste(sprintf("%s-vs-%s=%.4f",
                                combos[1,], combos[2,], smds), collapse = "; "))
}

# Binary/categorical SMD via Cohen's h
smd_binary <- function(var) {
  f <- as.formula(paste0("~", var))
  props <- svyby(f, ~bp_phenotype, pd, svymean, na.rm = TRUE)
  # take the TRUE / second column proportion
  pcol <- grep("TRUE|Yes|Male|Current", names(props), value = TRUE)[1]
  if (is.na(pcol)) pcol <- names(props)[2]
  pp <- setNames(props[[pcol]], props$bp_phenotype)
  combos <- combn(pheno_levels, 2)
  hs <- apply(combos, 2, function(pr) {
    a <- pr[1]; b <- pr[2]
    abs(2*asin(sqrt(pmin(pmax(pp[a],0),1))) - 2*asin(sqrt(pmin(pmax(pp[b],0),1))))
  })
  tibble(variable = var, type = "binary",
         max_smd = max(hs, na.rm = TRUE),
         detail = paste(sprintf("%s-vs-%s=%.4f",
                                combos[1,], combos[2,], hs), collapse = "; "))
}

smd_cont_vars <- c("age","bmi","egfr","hba1c","uacr","total_chol","hdl_chol",
                   "diet_sodium","antihtn_drug_classes","prevent_score")
smd_bin_vars  <- c("diabetes_dx","uninsured")
smd_tbl <- bind_rows(
  map_dfr(smd_cont_vars, ~ tryCatch(smd_continuous(.x), error = function(e) NULL)),
  map_dfr(smd_bin_vars,  ~ tryCatch(smd_binary(.x),     error = function(e) NULL))
) %>% mutate(meaningful_imbalance = max_smd > 0.1)
write_tbl(smd_tbl, "analysis2_smd_by_phenotype")

# =====================================================================
# ANALYSIS 3 — Phenotype predictors (multinomial via 2 binary Poisson)
# =====================================================================
log_step("Analysis 3: Phenotype predictors (two-model multinomial)")

covariate_rhs <- paste(
  "age10 + sex_f + race_eth_f + pir + education_f + uninsured +",
  "bmi5 + smoking_f + alcohol_f + egfr_cat_f + uacr_cat_f + hba1c +",
  "ckm_stage_f"
)

# Model A: Uncontrolled-Treated vs Controlled-Treated
dA <- subset(design, bp_phenotype %in% c("Controlled-Treated","Uncontrolled-Treated"))
dA <- update(dA, y = as.integer(bp_phenotype == "Uncontrolled-Treated"))
mA <- svyglm(as.formula(paste("y ~", covariate_rhs)), design = dA,
             family = quasipoisson(link = "log"))
write_tbl(pr_table(mA), "analysis3_modelA_uncontrolledTx_vs_controlled")

# Model B: Uncontrolled-Untreated vs Controlled-Treated
dB <- subset(design, bp_phenotype %in% c("Controlled-Treated","Uncontrolled-Untreated"))
dB <- update(dB, y = as.integer(bp_phenotype == "Uncontrolled-Untreated"))
mB <- svyglm(as.formula(paste("y ~", covariate_rhs)), design = dB,
             family = quasipoisson(link = "log"))
write_tbl(pr_table(mB), "analysis3_modelB_uncontrolledUntx_vs_controlled")

# =====================================================================
# ANALYSIS 4 — Primary regression: predictors of BP non-attainment
# Three sequential models, modified Poisson -> PRs
# =====================================================================
log_step("Analysis 4: Sequential PR models for BP non-attainment")

m1_f <- bp_nonattained ~ age10 + sex_f + race_eth_f
m2_f <- update(m1_f, . ~ . + pir + education_f + uninsured)
m3_f <- update(m2_f, . ~ . + bmi5 + smoking_f + alcohol_f +
                 egfr_cat_f + uacr_cat_f + hba1c + antihtn_drug_classes +
                 ckm_stage_f + diet_sodium_500 +
                 cotinine_log + pa_met_min_wk_1000)

# Cycle L lacks drug-class -> exclude L from any model containing
# antihtn_drug_classes (Model 3). Models 1-2 keep all cycles.
design_m3 <- subset(design_reg, nhanes_cycle != "L" & !is.na(antihtn_drug_classes))

m1 <- svyglm(m1_f, design = design_reg, family = quasipoisson(link = "log"))
m2 <- svyglm(m2_f, design = design_reg, family = quasipoisson(link = "log"))
m3 <- svyglm(m3_f, design = design_m3,  family = quasipoisson(link = "log"))

write_tbl(pr_table(m1) %>% mutate(model = "M1_demographics"),       "analysis4_model1_demographics")
write_tbl(pr_table(m2) %>% mutate(model = "M2_plus_SDOH"),          "analysis4_model2_sdoh")
write_tbl(pr_table(m3) %>% mutate(model = "M3_full_clinical"),      "analysis4_model3_full")

# Multicollinearity: unweighted OLS VIF on the full model formula
vif_dat <- design_m3$variables
vif_lm <- lm(update(m3_f, as.numeric(bp_nonattained) ~ .), data = vif_dat)
vif_out <- tryCatch({
  v <- car::vif(vif_lm)
  if (is.matrix(v)) {
    tibble(term = rownames(v), GVIF = v[,1], Df = v[,2],
           GVIF_adj = v[,3])  # GVIF^(1/(2Df))
  } else tibble(term = names(v), VIF = as.numeric(v))
}, error = function(e) tibble(note = paste("VIF failed:", conditionMessage(e))))
write_tbl(vif_out, "analysis4_vif")

# Supplementary: logistic ORs for comparison
m3_logit <- svyglm(m3_f, design = design_m3, family = quasibinomial(link = "logit"))
or_tbl <- pr_table(m3_logit) %>% rename(OR = PR)
write_tbl(or_tbl, "analysis4_supp_logistic_OR")

# =====================================================================
# ANALYSIS 5 — UACR × CKM Stage interaction (Bonferroni a = 0.017)
# =====================================================================
log_step("Analysis 5: UACR x CKM interaction")
int_f <- update(m3_f, . ~ . - uacr_cat_f - ckm_stage_f + uacr_cat_f * ckm_stage_f)
m_int <- svyglm(int_f, design = design_m3, family = quasipoisson(link = "log"))
p_int <- regTermTest(m_int, ~ uacr_cat_f:ckm_stage_f)$p
int_tbl <- tibble(
  test = "UACR x CKM stage",
  p_interaction = as.numeric(p_int),
  p_interaction_full = fmt_p(as.numeric(p_int)),
  bonferroni_alpha = 0.05/3,
  significant = as.numeric(p_int) < 0.05/3
)
write_tbl(int_tbl, "analysis5_uacr_ckm_interaction")

# =====================================================================
# ANALYSIS 6 — Temporal trend in mean SBP/DBP by cycle, per stage
# =====================================================================
log_step("Analysis 6: Temporal trends in mean BP")
design_trend <- update(design, cycle_year = dplyr::recode(nhanes_cycle,
                                                          G = 2012, H = 2014, I = 2016, J = 2018, L = 2022))

trend_rows <- map_dfr(c("sbp_mean","dbp_mean"), function(bpv) {
  map_dfr(levels(analytic$ckm_stage_f), function(st) {
    d_st <- subset(design_trend, ckm_stage_f == st)
    f <- as.formula(paste(bpv, "~ cycle_year"))
    mod <- tryCatch(svyglm(f, design = d_st, family = gaussian()),
                    error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    co <- summary(mod)$coefficients["cycle_year", ]
    tibble(outcome = bpv, ckm_stage = st,
           slope_mmHg_per_yr = co["Estimate"], se = co["Std. Error"],
           p_value = co[4], p_value_full = fmt_p(co[4]))
  })
})
write_tbl(trend_rows, "analysis6_temporal_bp_trends")

# Also cycle-stratified means (per CKM stage)
cyc_means <- svyby(~sbp_mean + dbp_mean, ~nhanes_cycle + ckm_stage_f,
                   design, svymean, na.rm = TRUE)
write_tbl(as_tibble(cyc_means), "analysis6b_cycle_stage_means")

# ---- Overall cycle-level means + overall trend (reference-format) ----------
# Matches the reference analysis9_temporal_bp_trend layout: one row per cycle
# with mean SBP/DBP, SE, 95% CI, and unweighted N; plus an overall linear trend
# (mmHg/yr) reported separately for SBP and DBP.
cyc_year_map <- c(G = 2012, H = 2014, I = 2016, J = 2018, L = 2022)
overall_means <- svyby(~sbp_mean + dbp_mean, ~nhanes_cycle, design,
                       svymean, na.rm = TRUE)
om_n <- analytic %>% count(nhanes_cycle, name = "n_unwt")
overall_tbl <- as_tibble(overall_means) %>%
  rename(sbp_mean = sbp_mean, dbp_mean = dbp_mean) %>%
  left_join(om_n, by = "nhanes_cycle") %>%
  mutate(
    year = cyc_year_map[nhanes_cycle],
    sbp_ci_lo = sbp_mean - 1.96 * se.sbp_mean, sbp_ci_hi = sbp_mean + 1.96 * se.sbp_mean,
    dbp_ci_lo = dbp_mean - 1.96 * se.dbp_mean, dbp_ci_hi = dbp_mean + 1.96 * se.dbp_mean
  ) %>%
  select(cycle = nhanes_cycle, year, n_unwt,
         sbp_mean, sbp_se = se.sbp_mean, sbp_ci_lo, sbp_ci_hi,
         dbp_mean, dbp_se = se.dbp_mean, dbp_ci_lo, dbp_ci_hi)
write_tbl(overall_tbl, "analysis6c_overall_cycle_means")

overall_trend <- map_dfr(c("sbp_mean","dbp_mean"), function(bpv) {
  mod <- svyglm(as.formula(paste(bpv, "~ cycle_year")),
                design = design_trend, family = gaussian())
  co <- summary(mod)$coefficients["cycle_year", ]
  tibble(outcome = bpv, slope_mmHg_per_yr = co["Estimate"], se = co["Std. Error"],
         p_value = co[4], p_value_full = fmt_p(co[4]),
         significant = co[4] < 0.05)
})
write_tbl(overall_trend, "analysis6c_overall_trend")

# ---- Per-phenotype cycle-level means + per-phenotype trend ------------------
# Restricted to the three treated/uncontrolled phenotypes (figure panels C/D).
pheno3 <- c("Controlled-Treated","Uncontrolled-Treated","Uncontrolled-Untreated")
design_phpheno <- subset(design_trend, bp_phenotype %in% pheno3)
pheno_cyc_means <- svyby(~sbp_mean + dbp_mean,
                         ~nhanes_cycle + bp_phenotype, design_phpheno,
                         svymean, na.rm = TRUE)
write_tbl(as_tibble(pheno_cyc_means), "analysis6d_cycle_phenotype_means")

pheno_trend <- map_dfr(c("sbp_mean","dbp_mean"), function(bpv) {
  map_dfr(pheno3, function(ph) {
    d_ph <- subset(design_phpheno, bp_phenotype == ph)
    mod <- tryCatch(svyglm(as.formula(paste(bpv, "~ cycle_year")),
                           design = d_ph, family = gaussian()),
                    error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    co <- summary(mod)$coefficients["cycle_year", ]
    tibble(outcome = bpv, phenotype = ph,
           slope_mmHg_per_yr = co["Estimate"], se = co["Std. Error"],
           p_value = co[4], p_value_full = fmt_p(co[4]),
           significant = co[4] < 0.05)
  })
})
write_tbl(pheno_trend, "analysis6d_phenotype_trend")

# =====================================================================
# ANALYSIS 7 — Fasting subgroup sensitivity
# =====================================================================
log_step("Analysis 7: Fasting subsample sensitivity")
if (!is.null(design_fast)) {
  design_fast2 <- update(design_fast,
                         glucose_fast_10 = glucose_fast / 10,
                         log_tg = log(pmax(triglycerides, 1)))
  f7 <- update(m3_f, . ~ . + glucose_fast_10 + log_tg)
  # exclude L (no drug classes) and require fasting weights
  dfa <- subset(design_fast2, nhanes_cycle != "L" & !is.na(antihtn_drug_classes))
  m7 <- tryCatch(svyglm(f7, design = dfa, family = quasipoisson(link = "log")),
                 error = function(e) NULL)
  if (!is.null(m7)) {
    write_tbl(pr_table(m7) %>% mutate(model = "fasting_sensitivity"),
              "analysis7_fasting_sensitivity")
  } else log_warn("  Analysis 7 model failed to converge; skipped.")
} else {
  log_warn("  Fasting design unavailable; Analysis 7 skipped.")
}

# =====================================================================
# ANALYSIS 8 — PREVENT statin sensitivity (3 scenarios x thresholds)
# =====================================================================
log_step("Analysis 8: PREVENT statin sensitivity")
cached_fn <- get("prevent_score_cached", envir = globalenv())
prevent_recalc <- function(statin_override) {
  d <- analytic
  if (statin_override == "primary") {
    # Primary uses actual statin status == the observed-data score -> cache hit.
    return(cached_fn(d))
  }
  # lower/upper override statin -> different inputs, bypass cache.
  st <- if (statin_override == "lower") rep(FALSE, nrow(d)) else rep(TRUE, nrow(d))
  cached_fn(d, statin_vec = st, use_cache = FALSE)
}
statin_scen <- tibble(scenario = c("lower","primary","upper")) %>%
  mutate(score = map(scenario, prevent_recalc))
thresholds <- c(7.5, 10)
statin_tbl <- map_dfr(seq_len(nrow(statin_scen)), function(i) {
  sc <- statin_scen$scenario[i]; scr <- statin_scen$score[[i]]
  dd <- update(design, .ps = scr)
  map_dfr(thresholds, function(th) {
    dd2 <- update(dd, exceed = as.integer(.ps >= th))
    m <- svymean(~exceed, dd2, na.rm = TRUE)
    tibble(scenario = sc, threshold = th,
           prop_exceed = as.numeric(m), se = as.numeric(SE(m)))
  })
})
write_tbl(statin_tbl, "analysis8_prevent_statin_sensitivity")

# =====================================================================
# ANALYSIS 9 — PREVENT HbA1c proxy
# =====================================================================
log_step("Analysis 9: PREVENT HbA1c proxy expansion")
proxy <- analytic %>% mutate(
  hba1c_proxy = dplyr::case_when(
    !is.na(hba1c) ~ hba1c,
    is.na(hba1c) & diabetes_dx %in% TRUE ~ 7.5,
    is.na(hba1c) & diabetes_dx %in% FALSE ~ 5.5,
    TRUE ~ NA_real_
  )
)
# Compare PREVENT-exceedance: complete-case (base model, cached) vs proxy set.
# The proxy only matters if HbA1c actually enters the equation, so the
# proxy-expanded score uses the PREVENT *hba1c* model with hba1c_proxy. The
# complete-case score is the cached base-model score (observed data).
cached_fn <- get("prevent_score_cached", envir = globalenv())
src <- get("prevent_risk", envir = globalenv())

proxy$prevent_score <- cached_fn(proxy)   # base-model, cached (complete-case)

# Proxy-expanded: hba1c model, scoring only rows where a proxy was applied
# (HbA1c originally missing); rows with observed HbA1c reuse the same hba1c-model
# score. We score the full set under the hba1c model once (no cache: different
# model + uses hba1c), which is the quantity of interest for this sensitivity.
proxy$prevent_proxy <- src(
  proxy$age, proxy$sex, proxy$sbp_mean, proxy$on_htn_meds,
  proxy$total_chol, proxy$hdl_chol,
  ifelse(is.na(proxy$on_statin), FALSE, proxy$on_statin),
  proxy$diabetes_dx, proxy$smoking_status == "Current",
  proxy$egfr, bmi = proxy$bmi, hba1c = proxy$hba1c_proxy,
  model = "hba1c"
)
d9 <- svydesign(ids=~sdmvpsu, strata=~sdmvstra, weights=~mec_wt_pooled,
                data = proxy, nest = TRUE)
cc <- svymean(~I(prevent_score >= 20), subset(d9, !is.na(prevent_score)), na.rm = TRUE)
px <- svymean(~I(prevent_proxy >= 20), subset(d9, !is.na(prevent_proxy)), na.rm = TRUE)
proxy_tbl <- tibble(
  set = c("complete_case","proxy_expanded"),
  prop_prevent_ge20 = c(as.numeric(cc)[1], as.numeric(px)[1]),
  se = c(as.numeric(SE(cc))[1], as.numeric(SE(px))[1])
)
write_tbl(proxy_tbl, "analysis9_prevent_hba1c_proxy")

# =====================================================================
# ANALYSIS 10 — Multiple imputation (MICE, PMM, M=10), parallelised
# PREVENT is re-scored on each imputed dataset (cannot be cached, since
# imputed bmi/egfr/lipids/hba1c differ per imputation). The 10 imputations
# are scored in PARALLEL via future.apply (multisession). Point estimate and
# SE are computed in a SINGLE pass per imputation (was two passes = 20 scorings;
# now 10). Live progress is written to logs/mice_progress.txt by the workers.
# =====================================================================
log_step("Analysis 10: MICE multiple imputation (M=10, PMM)")
N_CORES <- get("N_CORES", envir = globalenv())
M_IMP <- 10
mice_vars <- c("total_chol","hdl_chol","hba1c","egfr","bmi")
predictors <- c("age","sex","sbp_mean","diabetes_dx","smoking_status",
                "uacr","race_eth","education","pir")
mice_in <- analytic %>%
  select(all_of(unique(c(mice_vars, predictors,
                         "sdmvpsu","sdmvstra","mec_wt_pooled",
                         "bp_nonattained","on_htn_meds","ckm_stage_f")))) %>%
  mutate(across(where(is.character), as.factor))

# Parallel backend. Establish plan before calling futuremice.
future::plan(future::sequential)
future::plan(future::multisession, workers = min(N_CORES, M_IMP))
on.exit(future::plan(future::sequential), add = TRUE)

imp <- mice::futuremice(
  mice_in, m = M_IMP, method = "pmm",
  n.core = min(N_CORES, M_IMP), parallelseed = 2025
)

# Materialise the M completed datasets, trimmed to ONLY the columns scoring
# needs. This is the key performance fix: previously the whole `completed_list`
# was captured in the worker closure and exported in full to every worker
# (10x duplication, serialised per task). Now we (a) keep only needed columns,
# shrinking each frame, and (b) map future_lapply OVER the list so each worker
# receives exactly one dataset, not all ten.
score_cols <- c("age","sex","sbp_mean","on_htn_meds","total_chol","hdl_chol",
                "diabetes_dx","smoking_status","egfr","bmi",
                "sdmvpsu","sdmvstra","mec_wt_pooled")
completed_list <- lapply(seq_len(M_IMP), function(k) {
  ci <- mice::complete(imp, k)
  ci[, intersect(score_cols, names(ci)), drop = FALSE]
})

# Live-progress file the workers append to as each imputation finishes.
PROG_FILE <- file.path(PROJ_ROOT, "logs", "mice_progress.txt")
cat(sprintf("MICE PREVENT scoring progress (M=%d) — started %s\n",
            M_IMP, format(Sys.time())), file = PROG_FILE, append = FALSE)
log_step(sprintf("Scoring PREVENT on %d imputations in parallel (track: %s)",
                 M_IMP, basename(PROG_FILE)))
if (exists("notice")) notice(sprintf(
  "%d imputations scoring in parallel across %d workers. Live per-imputation lines stream to logs/%s (main console blocks until all finish, ~5-9 min).",
  M_IMP, min(N_CORES, M_IMP), basename(PROG_FILE)))

src <- get("prevent_risk", envir = globalenv())

# Score one imputed dataset (ci) -> estimate + SE. `k` is carried for logging.
score_one <- function(k, ci, src_fn, prog_path) {
  t0 <- Sys.time()
  ps <- src_fn(ci$age, ci$sex, ci$sbp_mean, ci$on_htn_meds %in% TRUE,
               ci$total_chol, ci$hdl_chol, FALSE,
               ci$diabetes_dx %in% TRUE, ci$smoking_status == "Current",
               ci$egfr, bmi = ci$bmi)
  di <- survey::svydesign(ids = ~sdmvpsu, strata = ~sdmvstra,
                          weights = ~mec_wt_pooled,
                          data = cbind(ci, ps = ps), nest = TRUE)
  m  <- survey::svymean(~I(ps >= 20), subset(di, !is.na(ps)), na.rm = TRUE)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("[%s] imputation %d/%d done in %.1fs (est=%.4f)\n",
              format(Sys.time(), "%H:%M:%S"), k, 10, dt, as.numeric(m)[1]),
      file = prog_path, append = TRUE)
  list(est = as.numeric(m)[1], se = as.numeric(survey::SE(m))[1])
}

# (Parallel backend established earlier, before futuremice)

# Map OVER the list (future_Map): each worker gets one frame + its index,
# never the whole list. src is a single small function (auto-exported once).
results <- future.apply::future_Map(
  function(k, ci) score_one(k, ci, src, PROG_FILE),
  seq_len(M_IMP), completed_list,
  future.seed = TRUE,
  future.packages = c("survey", "preventr")
)
future::plan(future::sequential)

# Echo the per-imputation timings to the console (workers wrote them to file
# during the parallel run; surface them now that the main process has them).
if (file.exists(PROG_FILE)) {
  cat("    --- MICE imputation timings ---\n")
  cat(paste0("    ", readLines(PROG_FILE)[-1]), sep = "\n"); cat("\n")
}
imp_estimates <- vapply(results, `[[`, numeric(1), "est")
imp_se        <- vapply(results, `[[`, numeric(1), "se")
log_ok(sprintf("All %d imputations scored. Pooling via Rubin's rules.", M_IMP))

# Rubin's rules
qbar <- mean(imp_estimates)
ubar <- mean(imp_se^2)
b    <- var(imp_estimates)
total_var <- ubar + (1 + 1/M_IMP) * b
df_rubin <- M_IMP - 1
tcrit <- qt(0.975, df_rubin)
mice_tbl <- tibble(
  metric = "PREVENT >= 20% prevalence (pooled)",
  pooled_estimate = qbar,
  total_SE = sqrt(total_var),
  ci_low = qbar - tcrit * sqrt(total_var),
  ci_high = qbar + tcrit * sqrt(total_var),
  within_var = ubar, between_var = b, df = df_rubin
)
write_tbl(mice_tbl, "analysis10_mice_pooled")

# ---- Exclusion flow summary -------------------------------------------------
flow <- readRDS(file.path(PROJ_ROOT, "data", "exclusion_flow.rds"))
write_tbl(tibble(stage = names(flow), n = as.integer(flow)), "analysis0_exclusion_flow")

log_ok("All analyses complete. Tables written to /output.")