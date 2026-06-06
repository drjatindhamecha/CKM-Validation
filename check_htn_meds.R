library(dplyr)
bpq_l <- readRDS("c:/Users/jatin/Downloads/CKM/CKM-Validation/data/BPQ_L.rds")
told_col="BPQ020"
med_col_l="BPQ101D"
med_col_gj=character(0)
old_htn <- dplyr::case_when(
  length(med_col_gj)>0 ~ bpq_l[[if(length(med_col_gj)) med_col_gj else told_col]]==1,
  length(med_col_l)>0 ~ bpq_l[[if(length(med_col_l)) med_col_l else told_col]]==1,
  TRUE~NA
)
new_htn <- if (length(med_col_gj)>0) {
  bpq_l[[med_col_gj]]==1
} else if (length(med_col_l)>0) {
  bpq_l[[med_col_l]]==1
} else {
  NA
}
cat("Cross-tabulation of Old vs New on_htn_meds logic in Cycle L:\n")
print(table(old=old_htn, new=new_htn, useNA="always"))
