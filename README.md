# CLIF National Pollution-Attributable ARF Geospatial Estimation Project

## CLIF VERSION  

2.1

---

## Objective

This project estimates the national burden of acute respiratory failure (ARF) using incomplete ICU network data from the CLIF consortium. Because CLIF captures ICU admissions from 44 hospitals across 10 health systems, observed ARF counts represent a subset of all ICU events nationally. We implement a Bayesian coverage-corrected geospatial modeling framework to:

1. Estimate county-level ICU admission coverage by CLIF  
2. Infer true county-level ARF incidence rates per capita  
3. Quantify geographic variation in ARF burden across the United States  
4. Estimate ARF cases and deaths attributable to long-term exposure to PM2.5 and NO2  

This framework treats ICU utilization and ARF incidence as latent processes and explicitly accounts for incomplete ICU coverage using county-of-residence linkage and all ICU admissions as an anchoring process.

Final outputs include county-level ARF incidence estimates with uncertainty, national burden estimates, and pollution-attributable ARF burden summaries.

---

## Required CLIF tables and fields

Please refer to the [CLIF data dictionary](https://clif-icu.com/data-dictionary), [CLIF Tools](https://clif-icu.com/tools), [ETL Guide](https://clif-icu.com/etl-guide), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for additional information.

The following CLIF tables are required:

### 1. **patient**
- `patient_id`
- `race_category`
- `ethnicity_category`
- `sex_category`
- `zip_code` or geographic identifier enabling county assignment  

**Rationale:** Required to assign county of residence and perform demographic poststratification using ACS distributions.

---

### 2. **hospitalization**
- `patient_id`
- `hospitalization_id`
- `admission_dttm`
- `discharge_dttm`
- `age_at_admission`

**Rationale:** Defines ICU admissions, enables annual aggregation, and supports age-based stratification.

---

### 3. **location**
- `hospitalization_id`
- `location_category`
- `location_start_dttm`
- `location_end_dttm`

**Rationale:** Identifies ICU encounters and determines ICU admission counts by county-year.

---

### 4. **respiratory_support**
- `hospitalization_id`
- `recorded_dttm`
- `device_category`
- `mode_category`
- `fio2_set`
- `peep_set`
- `resp_rate_set`
- `resp_rate_obs`

**Rationale:** Used to define physiologically grounded ARF criteria based on oxygen and ventilatory support.

---

### 5. **vitals**
- `hospitalization_id`
- `recorded_dttm`
- `vital_category`
- `vital_value`

Key categories:
- `resp_rate`
- `spo2`
- `map`
- `sbp`
- `dbp`

**Rationale:** Supports ARF phenotyping and severity characterization.

---

### 6. **labs**
- `hospitalization_id`
- `lab_result_dttm`
- `lab_category`
- `lab_value`

Key categories:
- `lactate`

**Rationale:** Used for severity adjustment and sensitivity analyses.

---

### 7. **discharge disposition / mortality**
- `hospitalization_id`
- in-hospital mortality indicator

**Rationale:** Required for ARF mortality modeling and pollution-attributable death estimation.

---

## External Data Requirements (Non-CLIF)

This project additionally requires:

- County-level annual population and demographic distributions (ACS)
- County-level annual average PM2.5 concentrations
- County-level annual average NO2 concentrations
- Optional: urbanicity, hospital density, distance to CLIF site (for coverage modeling)

---

## Cohort Identification

### Inclusion Criteria
- Adult ICU admissions occurring between January 1, 2018 and December 31, 2024
- Valid county-of-residence assignment
- ICU encounters identified using CLIF location data

### ARF Definition
ARF will be defined using physiologically grounded criteria incorporating:
- Invasive mechanical ventilation  
- Non-invasive ventilation  
- High-flow oxygen support  
- Elevated oxygen requirements and respiratory parameters  

Diagnosis codes alone will not define ARF.

### Exclusion Criteria
- Missing key demographic or geographic identifiers
- Non-ICU encounters
- Pediatric admissions (if applicable)

---

## Expected Results

The final project outputs will include:

1. Annual county-level ARF incidence rates per capita (posterior means and 95% credible intervals)
2. Annual county-level ICU coverage estimates by CLIF
3. National and regional ARF burden estimates
4. Maps of ARF incidence and uncertainty
5. Estimates of excess ARF cases attributable to PM2.5 and NO2
6. Estimates of pollution-attributable ARF mortality (secondary analysis)

All final outputs will be saved in: [output/final]

Intermediate modeling objects (INLA outputs, posterior samples) will be stored in: [output/models]


---

## Detailed Instructions for Running the Project

---

## 1. Update `config/config.json`

Follow instructions in the [config/README.md](config/README.md) file to specify:

- Data paths
- Study years
- Pollution exposure files
- Output directories
- Model options (e.g., SPDE mesh resolution, covariate inclusion)

**Note:** If using `01_run_cohort_id_app.R`, this step is optional.

---

## 2. Set Up the Project Environment

### R Environment (Primary)

Run: [code/00_renv_restore.R]


This will restore all required packages, including:
- INLA
- sf
- tidyverse
- data.table
- spdep
- terra

---


---

## 3. Run Code

The workflow proceeds in the following order:

1. **Cohort Identification**
   - Extract ICU admissions
   - Assign county of residence
   - Identify ARF cases

2. **County-Year Aggregation**
   - Compute total ICU admissions by county-year
   - Compute ARF admissions by county-year
   - Compute ARF mortality counts

3. **Merge External Data**
   - ACS population
   - PM2.5
   - NO2
   - Coverage predictors

4. **Modeling**
   - Fit ICU coverage model
   - Fit ARF incidence model
   - Fit ARF mortality model (optional)
   - Generate posterior summaries

5. **Attributable Burden Estimation**
   - Counterfactual pollution prediction
   - Excess ARF cases
   - Excess ARF deaths

Detailed instructions are provided in: [code/README.md]

---

## Example Repositories

- [CLIF Adult Sepsis Events](https://github.com/08wparker/CLIF_sepsis)
- [CLIF Eligibility for Mobilization](https://github.com/kaveriC/CLIF-eligibility-for-mobilization)
- [CLIF Variation in Ventilation](https://github.com/ingra107/clif_vent_variation)

---

## Summary

This repository implements a novel coverage-corrected geospatial modeling framework to extend incomplete ICU network data to national estimates of ARF burden. By jointly modeling ICU coverage, ARF incidence, and environmental exposures using Bayesian INLA-SPDE methods, this project produces county-level and national estimates of ARF and its attribution to long-term air pollution exposure.








