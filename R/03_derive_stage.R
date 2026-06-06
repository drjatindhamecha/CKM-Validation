# =============================================================================
# SEGMENT 3 — DERIVE VARIABLES & ASSIGN CKM STAGES
# Computes eGFR (2021 CKD-EPI race-free), UACR categories, PREVENT 10-yr CVD
# risk (with winsorisation + statin handling), BP attainment & phenotype,
# then applies the hierarchical CKM staging (4 -> 3 -> 2 -> 1 -> 0, higher
# stage preferred on overlap). Restricts the analytic sample to the study
# population (Stage 2+) per inclusion/exclusion rules.
# =============================================================================

stopifnot(exists("log_step"))
suppressPackageStartupMessages({ library(dplyr); library(preventr) })

merged_all <- get("merged_all", envir = globalenv())
df <- merged_all

# ---- 1. eGFR: 2021 CKD-EPI creatinine, race-free ---------------------------
# eGFR = 142 * min(Scr/k,1)^a * max(Scr/k,1)^-1.200 * 0.9938^Age * (1.012 if female)
# Note: the 2021 equation includes a 1.012 female multiplier. The spec's formula
# omitted it; the published Inker 2021 equation includes it. Included here for
# correctness; remove the female_mult line to match the spec's literal formula.
egfr_ckdepi <- function(scr, age, sex) {
  k <- ifelse(sex == "Female", 0.7, 0.9)
  a <- ifelse(sex == "Female", -0.241, -0.302)
  female_mult <- ifelse(sex == "Female", 1.012, 1.0)
  ratio <- scr / k
  142 * pmin(ratio, 1)^a * pmax(ratio, 1)^(-1.200) * (0.9938^age) * female_mult
}
df <- df %>% mutate(egfr = egfr_ckdepi(creatinine_s, age, sex))

# ---- 2. UACR category -------------------------------------------------------
df <- df %>% mutate(
  uacr_cat = dplyr::case_when(
    is.na(uacr)   ~ NA_character_,
    uacr < 30     ~ "A1",
    uacr < 300    ~ "A2",
    TRUE          ~ "A3"
  )
)

# ---- 3. Fasting flag (>= 8.5 hours) ----------------------------------------
df <- df %>% mutate(is_fasting = !is.na(fast_hours) & fast_hours >= 8.5)

# ---- 4. PREVENT 10-year total CVD risk (via preventr package) --------------
# Uses the official AHA PREVENT equations (Khan 2023) implemented in `preventr`
# rather than hand-coded coefficients. We winsorise inputs to PREVENT-valid
# ranges first (so ages >=80 -> 79, etc.) per the CKM staging protocol, then
# call estimate_risk() in vectorised data-frame mode and extract 10-year
# total_cvd risk as a percentage.
#
# NOTE on bounds: preventr internally caps SBP at 180 mmHg (AHA hypertensive-
# urgency rationale). We pre-winsorise to the published PREVENT ranges; for SBP
# we cap at 180 to match the package and avoid an out-of-range refusal. This is
# the staging-relevant choice used in the CKM literature (Lloyd-Jones/Ndumele).
winsorise <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# prevent_risk(): wrapper preserving the original call signature so the
# sensitivity analyses (statin scenarios, HbA1c proxy, MICE) keep working.
# `statin`, `diabetes`, `smoker` accept logical/0-1 vectors. Returns 10-yr
# total-CVD risk in PERCENT (NA for age < 30 or unusable inputs).
prevent_risk <- function(age, sex, sbp, on_bpmed, total_chol, hdl_chol,
                         on_statin, diabetes, smoker, egfr,
                         bmi = NULL, hba1c = NULL, uacr = NULL,
                         model = "base") {
  winsorise <- function(x, lo, hi) pmin(pmax(x, lo), hi)  # local: worker-safe
  n <- length(age)
  rec <- function(v) if (is.null(v)) rep(NA_real_, n) else v
  bmi <- rec(bmi); hba1c <- rec(hba1c); uacr <- rec(uacr)
  
  d <- data.frame(
    age     = winsorise(age, 30, 79),
    sex     = ifelse(sex == "Female", "female", "male"),
    sbp     = winsorise(sbp, 90, 180),
    bp_tx   = as.integer(ifelse(is.na(on_bpmed), FALSE, on_bpmed)),
    total_c = winsorise(total_chol, 130, 320),
    hdl_c   = winsorise(hdl_chol, 20, 100),
    statin  = as.integer(ifelse(is.na(on_statin), FALSE, on_statin)),
    dm      = as.integer(ifelse(is.na(diabetes), FALSE, diabetes)),
    smoking = as.integer(ifelse(is.na(smoker), FALSE, smoker)),
    egfr    = winsorise(egfr, 15, 150),
    bmi     = winsorise(bmi, 18.5, 40),
    stringsAsFactors = FALSE
  )
  if (!identical(model, "base")) {
    d$hba1c <- ifelse(is.na(hba1c), NA, winsorise(hba1c, 4.5, 15))
    d$uacr  <- uacr
  }
  
  req <- c("age","sbp","total_c","hdl_c","egfr","bmi")
  # Models that consume an optional predictor require it to be present, else
  # preventr would reject or silently downgrade the row.
  if (model %in% c("hba1c","full")) req <- c(req, "hba1c")
  if (model %in% c("uacr","full"))  req <- c(req, "uacr")
  ok  <- stats::complete.cases(d[, req]) & age >= 30
  out <- rep(NA_real_, n)
  if (!any(ok)) return(out)
  
  dd <- d[ok, , drop = FALSE]
  dd$preventr_id <- seq_len(nrow(dd))
  
  # add_to_dat = FALSE returns a clean frame keyed by `preventr_id` with a
  # numeric `total_cvd` column (verified against preventr 0.11 behaviour).
  res <- preventr::estimate_risk(
    use_dat = dd, model = model, time = "10yr",
    quiet = TRUE, collapse = TRUE, add_to_dat = FALSE, progress = FALSE
  )
  risk <- as.numeric(res$total_cvd) * 100
  pos  <- which(ok)
  out[pos[res$preventr_id]] <- risk
  out
}

# ---- SEQN-keyed PREVENT cache (with method-hash invalidation) --------------
# Observed-data PREVENT scores are a pure function of the inputs and the
# scoring method. We cache them keyed on SEQN in data/prevent_cache.rds so that
# Segment 3 (and Analyses 8/9 for unchanged rows) don't recompute the ~0.016s/
# row cost. A short hash of the scoring "method signature" is stored alongside;
# if the method changes (winsorisation bounds, model, wrapper edits), the hash
# changes and stale entries are ignored -- preventing silently-outdated scores.
PREVENT_CACHE <- file.path(PROJ_ROOT, "data", "prevent_cache.rds")
.prevent_method_hash <- function() {
  # Hash the wrapper body + key constants so any change invalidates the cache.
  sig <- paste(deparse(body(prevent_risk)), collapse = "\n")
  substr(digest::digest(sig, algo = "xxhash64"), 1, 16)
}

# prevent_score_cached(): compute PREVENT for a data frame `dat` that has SEQN
# plus the predictor columns, reusing cached scores where SEQN + method match.
# `statin_vec` lets callers override statin (Analysis 8); when overriding, the
# cache is bypassed for those rows (different inputs => different score).
prevent_score_cached <- function(dat, statin_vec = NULL, use_cache = TRUE) {
  meth <- tryCatch(.prevent_method_hash(), error = function(e) "nohash")
  n <- nrow(dat)
  scores <- rep(NA_real_, n)
  
  # Load cache (SEQN -> score) for the current method hash only.
  cache <- NULL
  if (use_cache && is.null(statin_vec) && file.exists(PREVENT_CACHE)) {
    cache <- tryCatch(readRDS(PREVENT_CACHE), error = function(e) NULL)
    if (!is.null(cache) && (is.null(cache$method) || cache$method != meth)) {
      cache <- NULL  # method changed -> discard stale cache
    }
  }
  
  hit <- rep(FALSE, n)
  if (!is.null(cache) && !is.null(cache$scores)) {
    idx <- match(dat$SEQN, cache$scores$SEQN)
    have <- !is.na(idx)
    scores[have] <- cache$scores$prevent_score[idx[have]]
    hit <- have & !is.na(scores)
  }
  
  need <- !hit
  if (any(need)) {
    sv <- if (is.null(statin_vec)) dat$on_statin else statin_vec
    fresh <- prevent_risk(
      dat$age[need], dat$sex[need], dat$sbp_mean[need], dat$on_htn_meds[need],
      dat$total_chol[need], dat$hdl_chol[need], sv[need],
      dat$diabetes_dx[need], dat$smoking_status[need] == "Current",
      dat$egfr[need], bmi = dat$bmi[need]
    )
    scores[need] <- fresh
    
    # Write freshly-computed observed-data scores back to cache (not overrides).
    if (use_cache && is.null(statin_vec)) {
      new_rows <- data.frame(SEQN = dat$SEQN[need], prevent_score = fresh)
      merged_scores <- if (!is.null(cache) && !is.null(cache$scores)) {
        rbind(cache$scores[!cache$scores$SEQN %in% new_rows$SEQN, , drop = FALSE],
              new_rows)
      } else new_rows
      saveRDS(list(method = meth, scores = merged_scores), PREVENT_CACHE)
    }
  }
  scores
}

# Compute PREVENT score for the cohort (base model), using the SEQN cache.
log_step("Computing PREVENT 10-yr total-CVD risk via preventr (cached by SEQN)")
if (exists("notice")) notice(sprintf(
  "Scoring PREVENT on ~%s participants in one pass — no live counter; expect ~5 min on first run (cached thereafter).",
  format(sum(df$age >= 30, na.rm = TRUE), big.mark = ",")))
df$prevent_score <- prevent_score_cached(df)
log_ok(sprintf("PREVENT scored for %d participants (non-NA)",
               sum(!is.na(df$prevent_score))))

# ---- 5. BP goal attainment & phenotype -------------------------------------
df <- df %>% mutate(
  bp_attained = !is.na(sbp_mean) & !is.na(dbp_mean) & sbp_mean < 130 & dbp_mean < 80,
  bp_nonattained = as.integer(!bp_attained),
  bp_phenotype = dplyr::case_when(
    is.na(sbp_mean) | is.na(dbp_mean) ~ NA_character_,
    on_htn_meds %in% TRUE  &  bp_attained ~ "Controlled-Treated",
    on_htn_meds %in% TRUE  & !bp_attained ~ "Uncontrolled-Treated",
    (!on_htn_meds %in% TRUE) & !bp_attained ~ "Uncontrolled-Untreated",
    (!on_htn_meds %in% TRUE) &  bp_attained ~ "Normotensive-Untreated",
    TRUE ~ NA_character_
  )
)

# aTRH per SAP: Uncontrolled-Treated & >=3 classes, OR Controlled-Treated & >=4
df <- df %>% mutate(
  atrh = dplyr::case_when(
    is.na(antihtn_drug_classes) ~ NA,
    bp_phenotype == "Uncontrolled-Treated" & antihtn_drug_classes >= 3 ~ TRUE,
    bp_phenotype == "Controlled-Treated"   & antihtn_drug_classes >= 4 ~ TRUE,
    TRUE ~ FALSE
  )
)

# ---- 6. CKM staging helpers -------------------------------------------------
asian <- function(re) re == "Non-Hispanic Asian"

bmi_elevated <- function(bmi, re) {
  thr <- ifelse(asian(re), 23, 25)
  !is.na(bmi) & bmi >= thr
}
waist_elevated <- function(waist, sex, re) {
  thr <- dplyr::case_when(
    asian(re) & sex == "Female" ~ 80,
    asian(re) & sex == "Male"   ~ 90,
    sex == "Female"             ~ 88,
    sex == "Male"               ~ 102,
    TRUE ~ NA_real_
  )
  !is.na(waist) & waist >= thr
}

# Stage 4: clinical CVD
stage4_crit <- function(d) d$cvd_history %in% TRUE

# Stage 3: high predicted risk OR very-high-risk CKD
stage3_crit <- function(d) {
  high_prevent <- !is.na(d$prevent_score) & d$prevent_score >= 20 &
    d$age >= 30 & d$age <= 79 & !(d$cvd_history %in% TRUE)
  ckd_vhigh <- (!is.na(d$egfr) & d$egfr < 30) |
    (!is.na(d$egfr) & d$egfr >= 30 & d$egfr < 60 & !is.na(d$uacr) & d$uacr >= 300) |
    (!is.na(d$uacr) & d$uacr >= 300)
  high_prevent | ckd_vhigh
}

# Stage 2: HTN | DM | high TG | low HDL | moderate-high CKD
stage2_crit <- function(d) {
  htn <- (!is.na(d$sbp_mean) & d$sbp_mean >= 130) |
    (!is.na(d$dbp_mean) & d$dbp_mean >= 80) |
    (d$on_htn_meds %in% TRUE)
  dm  <- (!is.na(d$hba1c) & d$hba1c >= 6.5) |
    (d$diabetes_dx %in% TRUE) |
    (!is.na(d$glucose_fast) & d$is_fasting & d$glucose_fast >= 126)
  hightg <- !is.na(d$triglycerides) & d$is_fasting & d$triglycerides >= 135
  lowhdl <- (!is.na(d$hdl_chol) & d$sex == "Male"   & d$hdl_chol < 40) |
    (!is.na(d$hdl_chol) & d$sex == "Female" & d$hdl_chol < 50)
  ckd_mod <- (!is.na(d$egfr) & d$egfr >= 60 & d$egfr < 90 & !is.na(d$uacr) & d$uacr >= 30) |
    (!is.na(d$egfr) & d$egfr >= 30 & d$egfr < 60 & (is.na(d$uacr) | d$uacr < 300)) |
    (!is.na(d$uacr) & d$uacr >= 30 & d$uacr < 300 & (!is.na(d$egfr) & d$egfr >= 60))
  htn | dm | hightg | lowhdl | ckd_mod
}

# Stage 1: adiposity (BMI/waist) OR prediabetes, without Stage 2+ criteria
stage1_crit <- function(d) {
  adipos <- bmi_elevated(d$bmi, d$race_eth) | waist_elevated(d$waist_cm, d$sex, d$race_eth)
  predm  <- (!is.na(d$hba1c) & d$hba1c >= 5.7 & d$hba1c < 6.5) |
    (!is.na(d$glucose_fast) & d$is_fasting & d$glucose_fast >= 100 & d$glucose_fast < 126)
  adipos | predm
}

# ---- 7. Assign stages hierarchically (higher preferred) --------------------
log_step("Assigning CKM stages")
df <- df %>% mutate(
  .s4 = stage4_crit(.),
  .s3 = stage3_crit(.),
  .s2 = stage2_crit(.),
  .s1 = stage1_crit(.),
  ckm_stage = dplyr::case_when(
    .s4 ~ 4L,
    .s3 ~ 3L,
    .s2 ~ 2L,
    .s1 ~ 1L,
    TRUE ~ 0L
  )
) %>% select(-.s4, -.s3, -.s2, -.s1)

stage_tab <- table(df$ckm_stage, useNA = "ifany")
log_ok(paste("CKM stage distribution (unweighted):",
             paste(names(stage_tab), stage_tab, sep = "=", collapse = ", ")))

# ---- 8. Apply inclusion / exclusion ----------------------------------------
log_step("Applying inclusion/exclusion criteria")
n0 <- nrow(df)

analytic <- df %>%
  filter(age >= 20) %>%                                   # adults >=20
  filter(!is.na(sbp_mean) & !is.na(dbp_mean)) %>%         # has BP
  filter(!(pregnant %in% TRUE)) %>%                       # not pregnant
  # Require lab data needed for staging / PREVENT
  filter(!is.na(creatinine_s), !is.na(uacr), !is.na(hba1c),
         !is.na(total_chol), !is.na(hdl_chol)) %>%
  filter(ckm_stage >= 2)                                  # study population

log_ok(sprintf("Sample after exclusions: %d (from %d merged)", nrow(analytic), n0))
flow <- c(
  merged = n0,
  age_ge20 = sum(df$age >= 20, na.rm = TRUE),
  has_bp = sum(df$age >= 20 & !is.na(df$sbp_mean) & !is.na(df$dbp_mean), na.rm = TRUE),
  final_stage2plus = nrow(analytic)
)
saveRDS(flow, file.path(PROJ_ROOT, "data", "exclusion_flow.rds"))

# ---- 9. Build modelling factors & scaled predictors ------------------------
analytic <- analytic %>% mutate(
  age10    = age / 10,
  bmi5     = bmi / 5,
  diet_sodium_500 = diet_sodium / 500,
  # Objective measures: cotinine (log1p due to right skew) and PA MET-min/wk
  cotinine_log = log1p(pmax(cotinine, 0)),
  pa_met_min_wk_1000 = pa_met_min_wk / 1000,
  sex_f    = factor(sex, levels = c("Female","Male")),
  race_eth_f = relevel(factor(race_eth), ref = "Non-Hispanic White"),
  education_f = relevel(factor(education,
                               levels = c("Less than high school","High school/GED","More than high school")),
                        ref = "More than high school"),
  uninsured = !(insured %in% TRUE),
  smoking_f = relevel(factor(smoking_status, levels = c("Never","Former","Current")), ref = "Never"),
  alcohol_cat = dplyr::case_when(
    is.na(drinks_yr) ~ NA_character_,
    !drinks_yr ~ "None",
    sex == "Female" & drinks_day <= 1 ~ "Moderate",
    sex == "Male"   & drinks_day <= 2 ~ "Moderate",
    drinks_yr & is.na(drinks_day) ~ "Moderate",
    TRUE ~ "Heavy"
  ),
  alcohol_f = relevel(factor(alcohol_cat, levels = c("None","Moderate","Heavy")), ref = "None"),
  egfr_cat = dplyr::case_when(
    egfr >= 60 ~ ">=60",
    egfr >= 30 ~ "30-59",
    TRUE ~ "<30"
  ),
  egfr_cat_f = relevel(factor(egfr_cat, levels = c(">=60","30-59","<30")), ref = ">=60"),
  uacr_cat_f = relevel(factor(uacr_cat, levels = c("A1","A2","A3")), ref = "A1"),
  ckm_stage_f = relevel(factor(ckm_stage, levels = c(2,3,4)), ref = "2"),
  nlr = ifelse(!is.na(neut_abs) & !is.na(lymph_abs) & lymph_abs > 0,
               neut_abs / lymph_abs, NA_real_),
  sii = ifelse(!is.na(platelets) & !is.na(neut_abs) & !is.na(lymph_abs) & lymph_abs > 0,
               platelets * neut_abs / lymph_abs, NA_real_)
)

saveRDS(analytic, file.path(PROJ_ROOT, "data", "analytic.rds"))
assign("analytic", analytic, envir = globalenv())
assign("df_full", df, envir = globalenv())
log_ok(sprintf("Saved data/analytic.rds (%d rows, %d cols)", nrow(analytic), ncol(analytic)))