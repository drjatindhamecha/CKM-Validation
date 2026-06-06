log_step <- function(x) cat(x, "\n")
log_ok <- function(x) cat(x, "\n")
library(dplyr)
df <- readRDS("c:/Users/jatin/Downloads/CKM/CKM-Validation/data/merged_all.rds")
cat("Merged: ", nrow(df), "\n")
df <- df %>% filter(age >= 20)
cat("Age >= 20: ", nrow(df), "\n")
df <- df %>% filter(!is.na(sbp_mean) & !is.na(dbp_mean))
cat("Has BP: ", nrow(df), "\n")
df <- df %>% filter(!(pregnant %in% TRUE))
cat("Not pregnant: ", nrow(df), "\n")
df <- df %>% filter(!is.na(creatinine_s), !is.na(uacr), !is.na(hba1c), !is.na(total_chol), !is.na(hdl_chol))
cat("Has labs: ", nrow(df), "\n")

# Need to run ckm_stage derivation from 03_derive_stage.R
merged_all <- df
PROJ_ROOT <- "c:/Users/jatin/Downloads/CKM/CKM-Validation"
source("c:/Users/jatin/Downloads/CKM/CKM-Validation/R/03_derive_stage.R")
