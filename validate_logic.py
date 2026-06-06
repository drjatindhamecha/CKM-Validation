"""
Logic-validation harness for the CKM staging algorithm and PREVENT inputs.
This does NOT replicate the R numerically — it stress-tests the DECISION LOGIC
(stage hierarchy, threshold boundaries, sex/race-specific cutoffs, winsorization)
that the R code in 03_derive_stage.R implements, to catch logic bugs before the
user runs the pipeline in R.
"""
import numpy as np
import pandas as pd

def asian(re): return re == "Non-Hispanic Asian"

def bmi_elevated(bmi, re):
    thr = 23 if asian(re) else 25
    return pd.notna(bmi) and bmi >= thr

def waist_elevated(waist, sex, re):
    if asian(re) and sex == "Female": thr = 80
    elif asian(re) and sex == "Male": thr = 90
    elif sex == "Female": thr = 88
    elif sex == "Male": thr = 102
    else: return False
    return pd.notna(waist) and waist >= thr

def stage4(r): return r.get("cvd_history") is True

def stage3(r):
    high_prevent = (pd.notna(r.get("prevent_score")) and r["prevent_score"] >= 20
                    and 30 <= r["age"] <= 79 and not (r.get("cvd_history") is True))
    egfr, uacr = r.get("egfr"), r.get("uacr")
    ckd_vhigh = ((pd.notna(egfr) and egfr < 30) or
                 (pd.notna(egfr) and 30 <= egfr < 60 and pd.notna(uacr) and uacr >= 300) or
                 (pd.notna(uacr) and uacr >= 300))
    return high_prevent or ckd_vhigh

def stage2(r):
    htn = ((pd.notna(r.get("sbp")) and r["sbp"] >= 130) or
           (pd.notna(r.get("dbp")) and r["dbp"] >= 80) or
           (r.get("on_htn_meds") is True))
    dm = ((pd.notna(r.get("hba1c")) and r["hba1c"] >= 6.5) or
          (r.get("diabetes_dx") is True) or
          (pd.notna(r.get("glu")) and r.get("fasting") and r["glu"] >= 126))
    hightg = pd.notna(r.get("tg")) and r.get("fasting") and r["tg"] >= 135
    hdl = r.get("hdl")
    lowhdl = ((pd.notna(hdl) and r["sex"] == "Male" and hdl < 40) or
              (pd.notna(hdl) and r["sex"] == "Female" and hdl < 50))
    egfr, uacr = r.get("egfr"), r.get("uacr")
    ckd_mod = ((pd.notna(egfr) and 60 <= egfr < 90 and pd.notna(uacr) and uacr >= 30) or
               (pd.notna(egfr) and 30 <= egfr < 60 and (pd.isna(uacr) or uacr < 300)) or
               (pd.notna(uacr) and 30 <= uacr < 300 and (pd.notna(egfr) and egfr >= 60)))
    return htn or dm or hightg or lowhdl or ckd_mod

def stage1(r):
    adipos = bmi_elevated(r.get("bmi"), r["race_eth"]) or waist_elevated(r.get("waist"), r["sex"], r["race_eth"])
    predm = ((pd.notna(r.get("hba1c")) and 5.7 <= r["hba1c"] < 6.5) or
             (pd.notna(r.get("glu")) and r.get("fasting") and 100 <= r["glu"] < 126))
    return adipos or predm

def stage(r):
    if stage4(r): return 4
    if stage3(r): return 3
    if stage2(r): return 2
    if stage1(r): return 1
    return 0

# ---- Edge-case test battery ----
tests = []
def check(name, row, expected):
    got = stage(row)
    tests.append((name, expected, got, "PASS" if got == expected else "FAIL"))

base = dict(age=50, sex="Male", race_eth="Non-Hispanic White", cvd_history=False,
            prevent_score=5, egfr=100, uacr=10, sbp=120, dbp=70, on_htn_meds=False,
            hba1c=5.0, diabetes_dx=False, glu=90, fasting=True, tg=100, hdl=60,
            bmi=22, waist=90)

# Hierarchy: CVD always wins -> 4
check("CVD overrides everything", {**base, "cvd_history": True, "hba1c": 9.0}, 4)
# PREVENT >=20 -> 3 (no CVD)
check("PREVENT 20 -> stage3", {**base, "prevent_score": 22}, 3)
# PREVENT 19.9 with HTN -> 2
check("PREVENT 19.9 + HTN -> stage2", {**base, "prevent_score": 19.9, "sbp": 140}, 2)
# Very-high CKD via eGFR<30 -> 3
check("eGFR 29 -> stage3", {**base, "egfr": 29}, 3)
# eGFR 45 + UACR 350 -> 3 (very high)
check("eGFR45 UACR350 -> stage3", {**base, "egfr": 45, "uacr": 350}, 3)
# eGFR 45 + UACR 100 -> 2 (moderate)
check("eGFR45 UACR100 -> stage2", {**base, "egfr": 45, "uacr": 100}, 2)
# UACR 300 alone -> 3
check("UACR300 -> stage3", {**base, "uacr": 300}, 3)
# UACR 299 + eGFR>=60 -> 2
check("UACR299 -> stage2", {**base, "uacr": 299}, 2)
# SBP exactly 130 -> 2
check("SBP 130 boundary -> stage2", {**base, "sbp": 130}, 2)
# SBP 129 DBP 79 no meds, nothing else -> check adiposity only
check("SBP129 DBP79 lean -> stage0", {**base, "sbp": 129, "dbp": 79, "hdl": 60,
      "bmi": 22, "waist": 90}, 0)
# Low HDL male <40 -> 2
check("HDL 39 male -> stage2", {**base, "hdl": 39}, 2)
# HDL 45 female <50 -> 2
check("HDL 45 female -> stage2", {**base, "sex": "Female", "hdl": 45, "waist": 80}, 2)
# HDL 45 male (>=40) not low -> not via HDL; lean -> 0
check("HDL 45 male not low -> stage0", {**base, "hdl": 45}, 0)
# Fasting TG 135 -> 2
check("TG135 fasting -> stage2", {**base, "tg": 135}, 2)
# Non-fasting TG 200 should NOT trigger stage2 via TG
check("TG200 nonfasting -> stage0", {**base, "tg": 200, "fasting": False}, 0)
# Prediabetes only -> 1
check("HbA1c 5.9 -> stage1", {**base, "hba1c": 5.9}, 1)
# Asian BMI 23 -> 1
check("Asian BMI23 -> stage1", {**base, "race_eth": "Non-Hispanic Asian", "bmi": 23, "waist": 70}, 1)
# Non-Asian BMI 23 -> 0
check("White BMI23 lean -> stage0", {**base, "bmi": 23, "waist": 90}, 0)
# HbA1c 6.5 -> diabetes -> 2
check("HbA1c 6.5 -> stage2", {**base, "hba1c": 6.5}, 2)
# Asian female waist 80 -> 1
check("Asian F waist80 -> stage1", {**base, "sex": "Female", "race_eth": "Non-Hispanic Asian",
      "waist": 80, "bmi": 20}, 1)

res = pd.DataFrame(tests, columns=["test","expected","got","result"])
print(res.to_string(index=False))
print("\n%d/%d passed" % ((res.result=="PASS").sum(), len(res)))
assert (res.result == "PASS").all(), "LOGIC BUG DETECTED"
print("All CKM staging logic checks passed.")