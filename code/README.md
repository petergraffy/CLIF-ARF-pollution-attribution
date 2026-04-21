## Code directory

This directory contains scripts for the project workflow. The general workflow consists of two main steps: cohort identification, and linkage/analysis.

### General Workflow

0. First, initialize your R environment using `00_renv_restore`.

1. Run the `01_cohort_identification.R` script. 
   This script should:
   - Apply inclusion and exclusion criteria
   - Select required fields from each table
   - Produce a CONSORT-style diagram for patient selection
   - Produce masked aggregated exports for:
     - county-year burden modeling
     - county-year age-sex post-stratification
     - pooled year age-sex race-ethnicity post-stratification
   - Automatically discover deterministic masking fragments in `keys/` if present

**Once all code has been run, upload the entire output/final folder as-is to Box.**

