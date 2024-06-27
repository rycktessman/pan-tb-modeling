# pan-tb-modeling
Code used to conduct modeling and cost analysis of novel pan-TB treatment regimens. 

The scripts in the "main" folder include those used in the health impact and economic analysis of pan-TB regimens (Projected Effects of a Pan-Tuberculosis Treatment Regimen: Modeling to Unpack the Health and Economic Benefits)

Order to run these scripts:

1. Run model (run_model.R) 
2. Estimate incidence from model output (estimate_incidence.R). 
3. Update model output based on changing cohort sizes from incidence estimates (run_final_calcs.R). 
4. Calculate price thresholds (calc_thresholds.R). 
5. Calculate means and confidence intervals of raw output (summarize_output.R)

All scripts are set up to run on a high-performance computing cluster. 

The scripts in the "inputs_for_resistance_analysis" folder contain code used to generate inputs to the resistance-focused analysis of pan-TB regimens (Exploring the potential of a pan-tuberculosis treatment regimen to drive the emergence of novel resistance)

Order to run these scripts:

1. Run model (run_model.R) - functionally the same as "run_model.R" in the "main" folder, but includes minor differences (scenario names, fewer sensitivity analyses, etc.)
2. Generate inputs for resistance-acquisition analysis (gen_inputs_for_resistance_analysis.R)

