# CKM × Blood Pressure Goal Attainment — Replication Pipeline

Replicates: *BP Goal Attainment Under the 2025 ACC/AHA Guideline Across CKM
Syndrome Stages 2–4* using pooled NHANES G,H,I,J (2011–2018) + L (2021–2023).

## How to run
```r
# From RStudio or terminal, in the project root:
Rscript run_study.R
```
The master script sources four segments in order, logging each step to
`logs/run_log.txt` and **halting on the first error**. Outputs land in `output/`.

Requires internet access for `nhanesA` to reach the CDC servers. The sandbox
this was authored in had no R and a restricted network, so the code was
**statically reviewed and its staging/PREVENT decision logic was validated via
a separate harness (`validate_logic.py`, 20/20 edge cases passing)** — but it
has **not been executed end-to-end against live NHANES**. Run it once in your
environment and skim the log; the per-segment design makes any failure point
obvious.

## Segment structure
1. `01_setup_packages.R` — initialises `renv`, installs packages to the
   project library, writes `renv.lock` (R + package versions) and
   `logs/environment.txt` (incl. RStudio version if detectable).
2. `02_load_merge.R` — pulls every required table per cycle, harmonises across
   cycles, merges to one table on `SEQN`. Caches raw pulls as `data/*.rds`.
3. `03_derive_stage.R` — eGFR (2021 CKD-EPI race-free), UACR categories,
   PREVENT, BP phenotype/attainment, CKM staging, inclusion/exclusion, model
   factors.
4. `04_analyses.R` — survey designs + Analyses 1–10.

## Decisions made
- **Data acquisition:** `nhanesA::nhanes()` per cycle.
- **PREVENT:** computed via the **`preventr`** package (official AHA Khan 2023
  equations), base model, 10-year `total_cvd`. Replaces the earlier hand-coded
  coefficients. Inputs winsorised to PREVENT-valid ranges (SBP capped at 180 to
  match the package). Base model is forced (`model = "base"`) so staging cannot
  silently switch to the hba1c/uacr/full variants.
- **CKM Stage 3 threshold:** 10-year PREVENT ≥ **20%** (confirmed against the
  CKM literature; the ≥7.5% figure is the 2025 AHA/ACC *treatment* threshold,
  not a staging boundary).
- **aTRH definition:** SAP version — Uncontrolled-Treated & ≥3 classes, OR
  Controlled-Treated & ≥4 classes.
- **Cycle L drug classes:** `NA`; L excluded from any model/table using
  `antihtn_drug_classes` (Model 3, fasting sensitivity, drug-class distribution).
- **Objective measures added to Model 3:** serum **cotinine** (`LBXCOT`,
  log1p-transformed) and **physical activity** MET-min/week from PAQ
  (per-1000 units). PAQ is questionnaire-based (no accelerometer data for these
  cycles); cotinine is a true biomarker.
- **LE8 removed** by your decision (NHANES proxy components too weak); not
  reintroduced.
- **Bonferroni family = 3 tests** at α = 0.0167: UACR×CKM, race×PIR,
  race×PREVENT.

## Segment structure
1. `01_setup_packages.R` — renv + packages (now includes `preventr`).
2. `02_load_merge.R` — NHANES pull + harmonise + merge (now includes COT, PAQ).
3. `03_derive_stage.R` — eGFR, UACR cats, PREVENT (via preventr), phenotypes,
   CKM staging, inclusion/exclusion, model factors (now incl. cotinine_log,
   pa_met_min_wk_1000).
4. `04_analyses.R` — survey designs + Analyses 1–10.
5. `05_readded_analyses.R` — race-stratified Model 3, race×PIR & race×PREVENT
   interactions, UACR×CKM cell prevalence + stratified PRs, drug-class
   distribution + aTRH prevalence.

## New output files (Segment 5)
`analysis5A_race_stratified_regression`, `analysis5B_race_pir_interaction`,
`analysis5C_race_prevent_interaction`, `analysis5D_uacr_ckm_prevalence`,
`analysis5D_uacr_ckm_stratified_pr`, `analysis5E_drug_class_distribution`,
`analysis5E_atrh_prevalence`. Segment 5 also rewrites
`analysis5_uacr_ckm_interaction.csv` to apply the shared α = 0.0167 verdict.

## Points you should still verify
1. **PA MET-minutes mapping.** PAQ variable names (PAQ610/PAD615 etc.) are
   stable across G–L but confirm against the codebook for your `nhanesA`
   version; the loader falls back to NA where fields are absent.
2. **Analysis 9 (HbA1c proxy)** now interacts with PREVENT only through the
   diabetes flag, since the base model doesn't take HbA1c. Still valid, but
   interpret accordingly.
3. **eGFR female ×1.012 multiplier** is still included (per published 2021
   equation); remove the `female_mult` line in `03_derive_stage.R` to match a
   spec that omits it.

## Original-implementation note
The earlier hand-coded PREVENT produced scores ~3–6 percentage points high,
which inflated Stage 3 membership and was the main driver of the Stage 3
divergence found in `pipeline_reconciliation.csv`. Switching to `preventr`
should bring Stage 3 attainment into line with the reference pipeline; re-run
and compare `analysis1_primary_bp_attainment_by_stage.csv` to confirm.

## Output files
All tables are CSV with full untruncated p-values (`p_value_full`), SEs, and
SMDs where applicable, plus NCHS suppression flags. Filenames map to analyses
`analysis1_*` … `analysis10_*`, plus `analysis0_exclusion_flow`.
