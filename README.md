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

1. **patient**
   - `patient_id`
   - `birth_date`
   - `sex_category`
   - `race_category`
   - `ethnicity_category`
   - (optional) `preferred_language`

   **Why:** demographic completeness checks; age computation fallback.

2. **hospitalization**
   - `patient_id`
   - `hospitalization_id`
   - `admission_dttm`
   - `discharge_dttm`
   - `age_at_admission`
   - residence geography fields (site dependent; any of the following are sufficient):
     - `county_code`
     - `census_tract`
     - `zipcode_five_digit` / `zipcode_nine_digit`

   **Why:** anchors the unit of analysis and enables county-of-residence assignment.

3. **adt** (or equivalent location history table)
   - `hospitalization_id`
   - `in_dttm`
   - `out_dttm`
   - `location_category` (must support ICU identification)

   **Why:** identifies ICU entry time (`first_icu_in`), ICU LOS, and ICU cohort membership.

4. **respiratory_support**
   - `hospitalization_id`
   - `recorded_dttm`
   - `device_category`
   - `mode_category`
   - `fio2_set`

   **Why:** supports FiO₂ pairing and ARF clinical definition logic.

5. **vitals**
   - `hospitalization_id`
   - `recorded_dttm`
   - `vital_category` (needs `spo2`)
   - `vital_value`

   **Why:** continuous SpO₂ density rule; hypoxemia on room air.

6. **labs**
   - `hospitalization_id`
   - `lab_result_dttm`
   - `lab_category` (needs `po2_arterial`, `pco2_arterial`, `ph_arterial`)
   - numeric value field (site dependent; e.g., `lab_value_numeric`)

   **Why:** PaO₂/FiO₂ ratio; hypercapnia pairing (pCO₂ + pH).

7. **hospital_diagnosis**
   - `hospitalization_id`
   - `diagnosis_code`

   **Why:** perioperative control cohort identification (J95.82–J95.84).

## External Data Requirements (Non-CLIF)

This project additionally requires:

- County-level annual population and demographic distributions (ACS)
- County-level annual average PM2.5 concentrations
- County-level annual average NO2 concentrations
- Optional: urbanicity, hospital density, distance to CLIF site (for coverage modeling)

---

## Cohort Identification

**Inclusion criteria**:

1.  Adult patients (≥18 years) admitted to ICU between 2018–2024.
2.  At least one of the following criteria of acute respiratory failure
    is met:
    -   Acute hypoxemic respiratory failure (any one of the following)
        -   SpO2 less than 90% on room air
        -   PaO2 of 60 mm Hg or less on room air
        -   PaO2–FiO2 ratio of 300 or less (on any amount of FiO2)
    -   Acute hypercapnic respiratory failure (both of the following)
        -   PaCO2 of 45 mm Hg or more AND
        -   Arterial pH \< 7.35
3.  Available ABG and/or continuous pulse oximetry data within ±24h of
    ICU admission.
4.  Residential census tract and county code for environmental data linkage.

***-\> Note that mixed hypoxic and hypercapnic respiratory failure is
common and should be accounted for.***

***-\> Also note that for SpO2 and PaO2, these numbers are directly
affected by supplemental oxygen (FiO2) via whatever delivery mechanism.
The definitions for ARF we choose for these values will be on room air
(21% FiO2). P/F ratio can define ARF even on supplemental oxygen.***

**Exclusion criteria** 
- Missing key demographic data (age, sex, race). 
- Hospitalizations \<24 hours in ICU. 
- Repeat ICU stays within same hospitalization (only first considered for primary analysis).

**Bibliography for definitions of ARF**

1\. Lagina, M. & Valley, T. S. Diagnosis and Management of Acute
Respiratory Failure. Critical Care Clinics 40, 235–253 (2024).

2\. Baldomero, A. K. et al. Effectiveness and Harms of High-Flow Nasal
Oxygen for Acute Respiratory Failure: An Evidence Report for a Clinical
Guideline From the American College of Physicians. Ann Intern Med 174,
952–966 (2021).

3\. RENOVATE Investigators and the BRICNet Authors et al. High-Flow
Nasal Oxygen vs Noninvasive Ventilation in Patients With Acute
Respiratory Failure: The RENOVATE Randomized Clinical Trial. JAMA 333,
875 (2025).

4\. Mirabile, V. S., Shebl, E., Sankari, A. & Burns, B. Respiratory
Failure in Adults. in StatPearls (StatPearls Publishing, Treasure Island
(FL), 2025).

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








