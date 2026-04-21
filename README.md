# CLIF National Pollution-Attributable ARF Geospatial Estimation Project

## CLIF VERSION  

2.1.0

---

## Objective

This project estimates the national burden of pollution-attributable acute respiratory failure (ARF) using incomplete ICU network data from the CLIF consortium. Because CLIF captures ICU admissions from 44 hospitals across 10 health systems, observed ARF counts represent a subset of all ICU events nationally. We implement a Bayesian coverage-corrected geospatial modeling framework to:

1. Estimate county-level ICU admission coverage by CLIF  
2. Infer true county-level ARF incidence rates per capita  
3. Quantify geographic variation in ARF burden across the United States  
4. Estimate ARF cases and deaths attributable to long-term exposure to PM2.5 and NO2  

This framework treats ICU utilization and ARF incidence as latent processes and explicitly accounts for incomplete ICU coverage using county-of-residence linkage and all ICU admissions as an anchoring process.

Final outputs from this work include county-level ARF incidence estimates with uncertainty, national burden estimates, and pollution-attributable ARF burden summaries.

---

## FOR CLIF RUN: Required CLIF tables and fields

The following CLIF tables are required:

1. **patient**
   - `patient_id`
   - `birth_date`
   - `sex_category`
   - `race_category`
   - `ethnicity_category`
   - (optional) `preferred_language`

2. **hospitalization**
   - `patient_id`
   - `hospitalization_id`
   - `admission_dttm`
   - `discharge_dttm`
   - `age_at_admission`
   - residence geography fields:
     - `county_code`

3. **adt** 
   - `hospitalization_id`
   - `in_dttm`
   - `out_dttm`
   - `location_category`

4. **respiratory_support**
   - `hospitalization_id`
   - `recorded_dttm`
   - `device_category`
   - `mode_category`
   - `fio2_set`

5. **vitals**
   - `hospitalization_id`
   - `recorded_dttm`
   - `vital_category`
   - `vital_value`

6. **labs**
   - `hospitalization_id`
   - `lab_result_dttm`
   - `lab_category`
   - `lab_value_numeric`

7. **hospital_diagnosis**
   - `hospitalization_id`
   - `diagnosis_code`

## NOT FOR CLIF SITES: External Data Requirements (Non-CLIF)

This project additionally requires:

- County-level annual population and demographic distributions (ACS)
- County-level annual average PM2.5 concentrations
- County-level annual average NO2 concentrations
- Optional: urbanicity, hospital density, distance to CLIF site (for coverage modeling)

But these tables are not necessary for the sites to have. They will be used during the analysis after the CLIF run.

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
4.  Residential county code for environmental data linkage.

**Exclusion criteria** 
- Missing key geographic data county code. 
- Hospitalizations \<24 hours in ICU. 

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

## Expected Project Results

The final project outputs will include:

1. Annual county-level ARF incidence rates per capita (posterior means and 95% credible intervals)
2. Annual county-level ICU coverage estimates by CLIF
3. National and regional ARF burden estimates
4. Maps of ARF incidence and uncertainty
5. Estimates of excess ARF cases attributable to PM2.5 and NO2
6. Estimates of pollution-attributable ARF mortality (secondary analysis)

---

## Detailed Instructions for CLIF Sites Running the Project

---

## 1. Update `config/config.json`

Follow instructions in the [config/README.md](config/README.md) file to specify:

- repo path
- tables path
- site name
- file type

---

## 2. Set Up the Project Environment

### R Environment

Run: [code/00_renv_restore.R]

This will restore all required packages.

---

## 3. Download Deterministic Masking Key Fragments

Before running the export script, download this site's deterministic masking key fragments from the external key-generation website and place them in the repo-level [keys](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARF-pollution-attribution/keys/README.md) folder.

The script expects exactly one file from each family:

- `key_county_year_Fragment_*.csv`
- `key_county_year_age_sex_Fragment_*.csv`
- `key_year_age_sex_race_ethnicity_Fragment_*.csv`

These fragments should not be renamed unless the filename still matches the same pattern.

---

## 4. Run Code

The workflow proceeds in the following order within `01_cohort_identification`:

1. **Cohort Identification**
   - Extract ICU admissions
   - Assign county of residence
   - Identify ARF cases and mortalities

2. **County-Year Aggregation**
   - Compute total ICU admissions by county-year
   - Compute ARF admissions by county-year
   - Compute ARF mortality counts
   - Apply deterministic masking using the key fragments in `keys/`
   - Compute masked post-stratification tables for age-sex and pooled demographic analyses
   
3. **CONSORT-style inclusion/exclusion table**

Detailed instructions are provided in: [code/README.md]

---

## Summary

This repository implements a novel coverage-corrected geospatial modeling framework to extend incomplete ICU network data to national estimates of ARF burden. By jointly modeling ICU coverage, ARF incidence, and environmental exposures using Bayesian INLA-SPDE methods, this project produces county-level and national estimates of ARF and its attribution to long-term air pollution exposure.






