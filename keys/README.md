# Keys Directory

This folder stores deterministic masking key fragments returned to a site from the external key-generation website.

## What Goes Here

Place the downloaded fragment CSV files directly in this folder.

Expected filename patterns:

- `key_county_year_Fragment_*.csv`
- `key_county_year_age_sex_Fragment_*.csv`
- `key_year_age_sex_race_ethnicity_Fragment_*.csv`

Examples:

- `key_county_year_Fragment_A.csv`
- `key_county_year_age_sex_Fragment_B.csv`
- `key_year_age_sex_race_ethnicity_Fragment_A.csv`

The site export code discovers these files automatically based on the filename pattern.

## Important Notes

- Do not rename fragment files unless the filename still matches the required pattern.
- Keep only the current site's fragment files in this folder when preparing a site release.
- The fragment files are externally created and should match the table dimensions exactly.
