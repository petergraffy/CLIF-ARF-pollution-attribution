## Code directory

This directory contains scripts for the project workflow. The general workflow consists of two main steps: cohort identification, and linkage/analysis.

### General Workflow

0. First, initialize your R environment using `00_renv_restore`.

1. Configure `config/config.json` for your site.

2. Download the deterministic masking key fragments for your site and place them in the repo-level `keys/` folder.
   Required fragment families:
   - `key_county_year_Fragment_*.csv`
   - `key_county_year_age_sex_Fragment_*.csv`
   - `key_year_age_sex_race_ethnicity_Fragment_*.csv`

3. Run the `01_cohort_identification.R` script. 
   This script should:
   - Apply inclusion and exclusion criteria
   - Select required fields from each table
   - Produce a CONSORT-style diagram for patient selection
   - Produce masked aggregated exports for:
     - county-year burden modeling
     - county-year age-sex post-stratification
     - pooled year age-sex race-ethnicity post-stratification
   - Automatically discover deterministic masking fragments in `keys/`
   - Error if a required fragment is missing or does not cover all observed rows

### Key workflow

The site script reads the key fragments directly from `keys/` using the masking repo naming convention. The masked exports written to `output/final/` are full scaffold tables, not just the observed rows from the site.

The three masked export families are:
- `site_county_year_*.csv`
- `site_county_year_age_sex_*.csv`
- `site_year_age_sex_race_ethnicity_*.csv`

**Once all code has been run, upload the entire output/final folder as-is to Box.**
