import os
import pandas as pd
import numpy as np

val_dir = r"c:\Users\jatin\Downloads\CKM\CKM-Validation\output"
att_dir = r"c:\Users\jatin\Downloads\CKM\CKM-BP-Attainment\output\tables"

def compare_csvs(val_file, att_file, key_cols, metric_cols):
    print(f"\n--- Comparing: {val_file} vs {att_file} ---")
    val_path = os.path.join(val_dir, val_file)
    att_path = os.path.join(att_dir, att_file)
    
    if not os.path.exists(val_path):
        print(f"Validation file missing: {val_file}")
        return
    if not os.path.exists(att_path):
        print(f"Attainment file missing: {att_file}")
        return
        
    df_val = pd.read_csv(val_path)
    df_att = pd.read_csv(att_path)
    
    # Try to align and compare
    try:
        if len(df_val) != len(df_att):
            print(f"Row count mismatch: Validation={len(df_val)}, Attainment={len(df_att)}")
            
        # Specific comparisons
        if "bp_attainment" in val_file:
            print("Validation Attainment:")
            print(df_val[['ckm_stage_f', 'prop_attained']].to_string(index=False))
            print("Attainment Attainment:")
            if 'pct' in df_att.columns:
                print(df_att[['group', 'pct']].to_string(index=False))
                
        elif "model3_full" in val_file:
            print("Validation Predictors (Race & UACR):")
            subset_val = df_val[df_val['term'].str.contains('race|uacr', case=False, na=False)]
            print(subset_val[['term', 'PR', 'CI_low', 'CI_high']].to_string(index=False))
            
            print("Attainment Predictors (Race & UACR):")
            if 'term' in df_att.columns:
                subset_att = df_att[df_att['term'].str.contains('race|uacr', case=False, na=False)]
                print(subset_att[['term', 'PR', 'CI_low', 'CI_high']].to_string(index=False))
                
        elif "trend" in val_file and "overall_trend" in val_file:
            print("Validation Trends:")
            print(df_val[['outcome', 'slope_mmHg_per_yr', 'p_value']].to_string(index=False))
            print("Attainment Trends:")
            if 'outcome' in df_att.columns:
                 print(df_att[['outcome', 'slope_mmHg_per_yr', 'p_value']].to_string(index=False))
            
    except Exception as e:
        print(f"Error comparing: {e}")

compare_csvs("analysis1_primary_bp_attainment_by_stage.csv", "analysis3_bp_attainment.csv", [], [])
compare_csvs("analysis4_model3_full.csv", "analysis4_model3_full.csv", [], [])
compare_csvs("analysis6c_overall_trend.csv", "analysis5_temporal_bp_trends.csv", [], [])
compare_csvs("analysis5F_pharmacotherapy_thresholds.csv", "analysis6_prevent_pharmacotherapy_thresholds.csv", [], [])
