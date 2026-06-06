# =============================================================================
# SEGMENT 2 — DOWNLOAD & MERGE NHANES VARIABLES
# Pulls each required table for cycles G,H,I,J,L via nhanesA, harmonises
# variable names/formats across cycles, and produces one merged table keyed
# on SEQN with a cycle indicator. Caches raw pulls to /data as .rds.
# =============================================================================

stopifnot(exists("log_step"))
suppressPackageStartupMessages({
  library(nhanesA); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

DATA_DIR <- file.path(PROJ_ROOT, "data")

# Cycle -> table suffix map. G..J are 2011-2018; L is 2021-2023.
CYCLES <- c(G = "G", H = "H", I = "I", J = "J", L = "L")
CYCLE_YEARS <- c(G = "2011-2012", H = "2013-2014", I = "2015-2016",
                 J = "2017-2018", L = "2021-2023")

# Tables we need by stem. BP differs: BPX (G-J) vs BPXO (L).
# Some stems are absent in some cycles; loader tolerates that.
TABLE_STEMS <- c("DEMO","BMX","BIOPRO","GHB","GLU","ALB_CR","TCHOL","HDL",
                 "TRIGLY","HSCRP","CBC","FASTQX","HIQ","BPQ","MCQ","DIQ",
                 "SMQ","ALQ","DR1TOT","DR2TOT","COT","PAQ")

# ---- Cached single-table fetch ---------------------------------------------
fetch_table <- function(stem, cyc) {
  suffix <- CYCLES[[cyc]]
  tbl_name <- paste0(stem, "_", suffix)
  cache <- file.path(DATA_DIR, paste0(tbl_name, ".rds"))
  if (file.exists(cache)) {
    return(readRDS(cache))
  }
  df <- tryCatch(
    suppressWarnings(nhanesA::nhanes(tbl_name, translated = FALSE)),
    error = function(e) NULL
  )
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    log_warn(sprintf("  table not available: %s", tbl_name))
    return(NULL)
  }
  names(df) <- toupper(names(df))
  saveRDS(df, cache)
  log_ok(sprintf("  fetched %s (%d rows, %d cols)", tbl_name, nrow(df), ncol(df)))
  df
}

# ---- BP fetch with BPX/BPXO harmonisation ----------------------------------
# Returns SEQN, sbp_mean, dbp_mean (mean of readings 2 & 3).
fetch_bp <- function(cyc) {
  if (cyc == "L") {
    bp <- fetch_table("BPXO", cyc)
    if (is.null(bp)) return(NULL)
    # Oscillometric: BPXOSY1..3 / BPXODI1..3
    sy <- intersect(c("BPXOSY2","BPXOSY3"), names(bp))
    di <- intersect(c("BPXODI2","BPXODI3"), names(bp))
    out <- bp %>%
      transmute(
        SEQN     = SEQN,
        sbp_mean = rowMeans(across(all_of(sy)), na.rm = TRUE),
        dbp_mean = rowMeans(across(all_of(di)), na.rm = TRUE)
      )
  } else {
    bp <- fetch_table("BPX", cyc)
    if (is.null(bp)) return(NULL)
    sy <- intersect(c("BPXSY2","BPXSY3"), names(bp))
    di <- intersect(c("BPXDI2","BPXDI3"), names(bp))
    out <- bp %>%
      transmute(
        SEQN     = SEQN,
        sbp_mean = rowMeans(across(all_of(sy)), na.rm = TRUE),
        dbp_mean = rowMeans(across(all_of(di)), na.rm = TRUE)
      )
  }
  # DBP of 0 is a measurement artefact (Korotkoff disappearance recorded as 0):
  out$dbp_mean[!is.na(out$dbp_mean) & out$dbp_mean == 0] <- NA_real_
  out$sbp_mean[is.nan(out$sbp_mean)] <- NA_real_
  out$dbp_mean[is.nan(out$dbp_mean)] <- NA_real_
  out
}

# ---- Antihypertensive drug-class count from RXQ_RX (G-J only) --------------
# Matches ATC classes C02/C03/C07/C08/C09 via drug-name keyword lists.
ANTIHTN_KEYWORDS <- list(
  C03 = c("hydrochlorothiazide","chlorthalidone","furosemide","torsemide",
          "bumetanide","spironolactone","eplerenone","amiloride","triamterene",
          "indapamide","metolazone","chlorothiazide"),
  C07 = c("metoprolol","atenolol","carvedilol","bisoprolol","propranolol",
          "nebivolol","labetalol","nadolol","sotalol","pindolol","acebutolol"),
  C08 = c("amlodipine","nifedipine","diltiazem","verapamil","felodipine",
          "nicardipine","isradipine","nisoldipine"),
  C09 = c("lisinopril","enalapril","ramipril","benazepril","quinapril",
          "fosinopril","perindopril","trandolapril","captopril","moexipril",
          "losartan","valsartan","irbesartan","candesartan","olmesartan",
          "telmisartan","azilsartan","eprosartan","aliskiren","sacubitril"),
  C02 = c("clonidine","hydralazine","minoxidil","methyldopa","guanfacine",
          "doxazosin","prazosin","terazosin","reserpine")
)
STATIN_KEYWORDS <- c("simvastatin","atorvastatin","rosuvastatin","pravastatin",
                     "lovastatin","fluvastatin","pitavastatin")

fetch_meds_GJ <- function(cyc) {
  rx <- fetch_table("RXQ_RX", cyc)
  if (is.null(rx)) return(NULL)
  drug_col <- intersect(c("RXDDRUG"), names(rx))
  if (!length(drug_col)) return(NULL)
  rx <- rx %>% mutate(drug = tolower(as.character(.data[[drug_col[1]]])))
  per_person <- rx %>%
    group_by(SEQN) %>%
    summarise(
      classes = {
        d <- drug
        hits <- vapply(names(ANTIHTN_KEYWORDS), function(cl) {
          any(str_detect(d, str_c(ANTIHTN_KEYWORDS[[cl]], collapse = "|")), na.rm = TRUE)
        }, logical(1))
        sum(hits)
      },
      on_statin = any(str_detect(drug, str_c(STATIN_KEYWORDS, collapse = "|")), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(SEQN,
              antihtn_drug_classes = as.integer(classes),
              on_statin = on_statin)
  per_person
}

# ---- Assemble one cycle -----------------------------------------------------
build_cycle <- function(cyc) {
  log_step(sprintf("Assembling cycle %s (%s)", cyc, CYCLE_YEARS[[cyc]]))
  
  demo <- fetch_table("DEMO", cyc)
  if (is.null(demo)) stop(sprintf("DEMO missing for cycle %s — cannot proceed", cyc))
  
  # DEMO core fields (names stable across G-L)
  d <- demo %>% transmute(
    SEQN, nhanes_cycle = cyc,
    age      = RIDAGEYR,
    sex      = dplyr::recode(as.integer(RIAGENDR), `1` = "Male", `2` = "Female"),
    race_eth = dplyr::recode(as.integer(RIDRETH3),
                             `1`="Mexican American", `2`="Other Hispanic",
                             `3`="Non-Hispanic White", `4`="Non-Hispanic Black",
                             `6`="Non-Hispanic Asian", `7`="Other/Multi"),
    educ_raw = if ("DMDEDUC2" %in% names(demo)) DMDEDUC2 else NA_real_,
    pir      = if ("INDFMPIR" %in% names(demo)) INDFMPIR else NA_real_,
    preg_demo= if ("RIDEXPRG" %in% names(demo)) RIDEXPRG else NA_real_,
    WTMEC2YR = WTMEC2YR,
    sdmvpsu  = SDMVPSU,
    sdmvstra = SDMVSTRA
  ) %>%
    mutate(education = dplyr::case_when(
      educ_raw %in% c(1,2) ~ "Less than high school",
      educ_raw == 3        ~ "High school/GED",
      educ_raw %in% c(4,5) ~ "More than high school",
      TRUE                 ~ NA_character_
    )) %>% select(-educ_raw)
  
  bp  <- fetch_bp(cyc)
  bmx <- fetch_table("BMX", cyc)
  bio <- fetch_table("BIOPRO", cyc)
  ghb <- fetch_table("GHB", cyc)
  glu <- fetch_table("GLU", cyc)
  alb <- fetch_table("ALB_CR", cyc)
  tc  <- fetch_table("TCHOL", cyc)
  hdl <- fetch_table("HDL", cyc)
  tg  <- fetch_table("TRIGLY", cyc)
  crp <- fetch_table("HSCRP", cyc)
  cbc <- fetch_table("CBC", cyc)
  fas <- fetch_table("FASTQX", cyc)
  hiq <- fetch_table("HIQ", cyc)
  bpq <- fetch_table("BPQ", cyc)
  mcq <- fetch_table("MCQ", cyc)
  diq <- fetch_table("DIQ", cyc)
  smq <- fetch_table("SMQ", cyc)
  alq <- fetch_table("ALQ", cyc)
  dr1 <- fetch_table("DR1TOT", cyc)
  dr2 <- fetch_table("DR2TOT", cyc)
  cot <- fetch_table("COT", cyc)
  paq <- fetch_table("PAQ", cyc)
  
  pick <- function(df, cols, newnames) {
    if (is.null(df)) return(NULL)
    have <- intersect(cols, names(df))
    if (!length(have)) return(NULL)
    out <- df[, c("SEQN", have), drop = FALSE]
    names(out) <- c("SEQN", newnames[match(have, cols)])
    out
  }
  
  bmx_s <- pick(bmx, c("BMXBMI","BMXWAIST"), c("bmi","waist_cm"))
  bio_s <- pick(bio, c("LBXSCR","LBXSUA"),   c("creatinine_s","uric_acid"))
  ghb_s <- pick(ghb, c("LBXGH"),             c("hba1c"))
  glu_s <- pick(glu, c("LBXGLU"),            c("glucose_fast"))
  tc_s  <- pick(tc,  c("LBXTC"),             c("total_chol"))
  hdl_s <- pick(hdl, c("LBDHDD"),            c("hdl_chol"))
  tg_s  <- pick(tg,  c("LBXTR"),             c("triglycerides"))
  crp_s <- pick(crp, c("LBXHSCRP"),          c("hscrp"))
  fas_s <- pick(fas, c("PHAFSTHR"),          c("fast_hours"))
  
  # UACR: urinary albumin (mg/L) / urinary creatinine (g/L) * 100
  alb_s <- NULL
  if (!is.null(alb)) {
    ualb_col <- intersect(c("URXUMA"), names(alb))     # albumin ug/mL? -> mg/L
    ucr_col  <- intersect(c("URXUCR"), names(alb))     # creatinine mg/dL
    if (length(ualb_col) && length(ucr_col)) {
      # Standard UACR (mg/g) = urinary albumin (mg/L) / urinary creatinine (g/L).
      # NHANES URXUMA is in ug/mL == mg/L; URXUCR is in mg/dL.
      # creatinine (g/L) = URXUCR(mg/dL) * 10 / 1000 = URXUCR * 0.01
      alb_s <- alb %>% transmute(
        SEQN,
        uacr = .data[[ualb_col]] / (.data[[ucr_col]] * 0.01)
      )
    }
  }
  
  # CBC for NLR / SII: neutrophils & lymphocytes (1000 cells/uL) + platelets
  cbc_s <- pick(cbc,
                c("LBXNEPCT","LBXLYPCT","LBDNENO","LBDLYMNO","LBXPLTSI"),
                c("neut_pct","lymph_pct","neut_abs","lymph_abs","platelets"))
  
  hiq_s <- NULL
  if (!is.null(hiq)) {
    icol <- intersect(c("HIQ011"), names(hiq))
    if (length(icol)) hiq_s <- hiq %>% transmute(SEQN, insured = .data[[icol]] == 1)
  }
  
  # BPQ — self-reported HTN dx + medication use. Harmonise across cycles.
  bpq_s <- NULL
  if (!is.null(bpq)) {
    told_col <- intersect(c("BPQ020"), names(bpq))                 # ever told high BP
    med_col_gj <- intersect(c("BPQ050A"), names(bpq))              # now taking BP meds (G-J)
    med_col_l  <- intersect(c("BPQ101D"), names(bpq))              # taking Rx for HBP (L)
    bpq_s <- bpq %>% transmute(
      SEQN,
      told_htn = if (length(told_col)) .data[[told_col]] == 1 else NA,
      on_htn_meds = dplyr::case_when(
        length(med_col_gj) > 0 ~ .data[[if (length(med_col_gj)) med_col_gj else told_col]] == 1,
        length(med_col_l)  > 0 ~ .data[[if (length(med_col_l)) med_col_l else told_col]] == 1,
        TRUE ~ NA
      )
    )
  }
  
  mcq_s <- NULL
  if (!is.null(mcq)) {
    cols <- intersect(c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"), names(mcq))
    if (length(cols)) {
      mcq_s <- mcq %>% transmute(
        SEQN,
        cvd_history = rowSums(across(all_of(cols), ~ .x == 1), na.rm = TRUE) > 0,
        mcq_hf  = if ("MCQ160B" %in% cols) MCQ160B == 1 else FALSE,
        mcq_chd = if ("MCQ160C" %in% cols) MCQ160C == 1 else FALSE,
        mcq_ang = if ("MCQ160D" %in% cols) MCQ160D == 1 else FALSE,
        mcq_mi  = if ("MCQ160E" %in% cols) MCQ160E == 1 else FALSE,
        mcq_str = if ("MCQ160F" %in% cols) MCQ160F == 1 else FALSE
      )
    }
  }
  
  diq_s <- NULL
  if (!is.null(diq)) {
    dcol <- intersect(c("DIQ010"), names(diq))
    if (length(dcol)) diq_s <- diq %>% transmute(SEQN, diabetes_dx = .data[[dcol]] == 1)
  }
  
  # Smoking: SMQ020 (ever 100 cigs) + SMQ040 (now smoke) -> never/former/current
  smq_s <- NULL
  if (!is.null(smq)) {
    if (all(c("SMQ020","SMQ040") %in% names(smq))) {
      smq_s <- smq %>% transmute(
        SEQN,
        smoking_status = dplyr::case_when(
          SMQ020 == 2 ~ "Never",
          SMQ020 == 1 & SMQ040 %in% c(1,2) ~ "Current",
          SMQ020 == 1 & SMQ040 == 3 ~ "Former",
          TRUE ~ NA_character_
        )
      )
    }
  }
  
  # Alcohol harmonisation: ALQ101 (G-I) vs ALQ111 (J-L); ALQ130 drinks/day
  alq_s <- NULL
  if (!is.null(alq)) {
    any_col <- intersect(c("ALQ101","ALQ111"), names(alq))
    qty_col <- intersect(c("ALQ130"), names(alq))
    if (length(any_col)) {
      alq_s <- alq %>% transmute(
        SEQN,
        drinks_yr = .data[[any_col[1]]] == 1,
        drinks_day = if (length(qty_col)) suppressWarnings(as.numeric(.data[[qty_col]])) else NA_real_
      )
    }
  }
  
  # Sodium: mean of DR1TSODI & DR2TSODI
  diet_s <- NULL
  d1 <- pick(dr1, c("DR1TSODI"), c("na1"))
  d2 <- pick(dr2, c("DR2TSODI"), c("na2"))
  if (!is.null(d1) || !is.null(d2)) {
    diet_s <- full_join(d1 %||% tibble(SEQN = integer()),
                        d2 %||% tibble(SEQN = integer()), by = "SEQN") %>%
      mutate(diet_sodium = rowMeans(across(any_of(c("na1","na2"))), na.rm = TRUE)) %>%
      mutate(diet_sodium = ifelse(is.nan(diet_sodium), NA_real_, diet_sodium)) %>%
      select(SEQN, diet_sodium)
  }
  
  meds_s <- if (cyc != "L") fetch_meds_GJ(cyc) else NULL
  
  # Cotinine (objective nicotine biomarker): serum cotinine ng/mL = LBXCOT
  cot_s <- NULL
  if (!is.null(cot)) {
    ccol <- intersect(c("LBXCOT"), names(cot))
    if (length(ccol)) cot_s <- cot %>% transmute(SEQN, cotinine = .data[[ccol]])
  }
  
  # Physical activity: total weekly MET-minutes from PAQ global items.
  # NHANES GPAQ-style: moderate/vigorous work + recreation minutes/week.
  # MET-min/week = sum(minutes * MET weight). Vigorous=8 MET, Moderate=4 MET.
  # Variables: PAQ605/610/615 (vigorous work days/min), PAQ620/625/630 (mod work),
  # PAQ650/655/660 (vigorous rec), PAQ665/670/675 (mod rec) vary by cycle; we
  # use the summary minute variables where present and fall back gracefully.
  paq_s <- NULL
  if (!is.null(paq)) {
    # Common continuous fields: PAD615 (vig work min/day), PAD630 (mod work),
    # PAD660 (vig rec), PAD675 (mod rec); day-count fields PAQ610/625/655/670.
    has <- function(v) v %in% names(paq)
    safe <- function(v) if (has(v)) suppressWarnings(as.numeric(paq[[v]])) else rep(NA_real_, nrow(paq))
    # minutes/week = days/week * minutes/day; guard implausible codes (>=7777)
    clamp <- function(x) ifelse(!is.na(x) & x >= 7777, NA_real_, x)
    vw <- clamp(safe("PAQ610")) * clamp(safe("PAD615"))   # vigorous work
    mw <- clamp(safe("PAQ625")) * clamp(safe("PAD630"))   # moderate work
    vr <- clamp(safe("PAQ655")) * clamp(safe("PAD660"))   # vigorous recreation
    mr <- clamp(safe("PAQ670")) * clamp(safe("PAD675"))   # moderate recreation
    met_min <- rowSums(cbind(8*vw, 4*mw, 8*vr, 4*mr), na.rm = TRUE)
    # If a person had no usable PA fields at all, set NA rather than 0
    all_na <- is.na(vw) & is.na(mw) & is.na(vr) & is.na(mr)
    met_min[all_na] <- NA_real_
    paq_s <- tibble(SEQN = paq$SEQN, pa_met_min_wk = met_min)
  }
  
  # Pregnancy from urine test if available (else DEMO RIDEXPRG)
  preg_s <- NULL
  if (!is.null(demo) && "preg_demo" %in% names(d)) {
    preg_s <- d %>% transmute(SEQN, pregnant = preg_demo == 1)
  }
  
  # Fasting subsample weights (for fasting analyses)
  fastwt_s <- NULL
  if (!is.null(tg)) {
    wcol <- intersect(c("WTSAF2YR"), names(tg))
    if (length(wcol)) fastwt_s <- tg %>% transmute(SEQN, WTSAF2YR = .data[[wcol]])
  }
  if (is.null(fastwt_s) && !is.null(glu)) {
    wcol <- intersect(c("WTSAF2YR"), names(glu))
    if (length(wcol)) fastwt_s <- glu %>% transmute(SEQN, WTSAF2YR = .data[[wcol]])
  }
  
  # ---- Reduce-join all components on SEQN -----------------------------------
  parts <- list(d, bp, bmx_s, bio_s, ghb_s, glu_s, alb_s, tc_s, hdl_s, tg_s,
                crp_s, cbc_s, fas_s, hiq_s, bpq_s, mcq_s, diq_s, smq_s, alq_s,
                diet_s, meds_s, cot_s, paq_s, preg_s, fastwt_s)
  parts <- Filter(Negate(is.null), parts)
  merged <- Reduce(function(a, b) left_join(a, b, by = "SEQN"), parts)
  
  # Pooled weights (5 cycles) per spec
  merged$mec_wt_pooled <- merged$WTMEC2YR / 5
  if ("WTSAF2YR" %in% names(merged)) merged$saf_wt_pooled <- merged$WTSAF2YR / 5
  
  # Cycle L lacks drug-class detail -> NA so it is excluded from those models
  if (cyc == "L") {
    merged$antihtn_drug_classes <- NA_integer_
    merged$on_statin <- NA
  }
  log_ok(sprintf("Cycle %s assembled: %d participants, %d columns",
                 cyc, nrow(merged), ncol(merged)))
  merged
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Build all cycles & row-bind -------------------------------------------
all_cycles <- lapply(names(CYCLES), build_cycle)

# Align columns (union) before binding
all_cols <- Reduce(union, lapply(all_cycles, names))
all_cycles <- lapply(all_cycles, function(df) {
  miss <- setdiff(all_cols, names(df))
  for (m in miss) df[[m]] <- NA
  df[, all_cols, drop = FALSE]
})
merged_all <- bind_rows(all_cycles)

log_ok(sprintf("MERGED ALL CYCLES: %d participants, %d variables",
               nrow(merged_all), ncol(merged_all)))

saveRDS(merged_all, file.path(DATA_DIR, "merged_all.rds"))
assign("merged_all", merged_all, envir = globalenv())
log_ok("Saved data/merged_all.rds")