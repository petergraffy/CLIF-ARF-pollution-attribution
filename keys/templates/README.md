## Key Templates

These files are meant to mirror the deterministic masking project more closely: the real artifacts should be complete scaffolds, not two-row examples.

Files in this folder:
- `generate_key_templates.py`: generates full scaffold CSVs
- `key_county_year_Fragment_TEMPLATE.csv`
- `key_county_year_age_sex_Fragment_TEMPLATE.csv`
- `key_year_age_sex_race_ethnicity_Fragment_TEMPLATE.csv`

General rules:
- Create one key CSV per site per output table.
- Keep the dimension names exactly as written in the scaffold headers.
- These templates now mirror the masking repo fragment style: join columns plus `offset`.
- Replace the placeholder `0` offsets with the actual site-specific offsets returned by the external system.

Expected fragment naming inside the project:
- `key_county_year_Fragment_<id>.csv`
- `key_county_year_age_sex_Fragment_<id>.csv`
- `key_year_age_sex_race_ethnicity_Fragment_<id>.csv`

Put the downloaded fragments in the repo-level [keys](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARF-pollution-attribution/keys/README.md) folder.

Allowed dimension values:
- `year`: `2018` through `2024`
- `age_band`: `18-39`, `40-64`, `65-74`, `75+`
- `sex`: `Female`, `Male`, `Other/Unknown`
- `race_group`: `AIAN`, `Asian`, `Black`, `NHPI`, `White`, `Other/Unknown`
- `ethnicity_group`: `Hispanic`, `Non-Hispanic`, `Other`

Note on counties:
- The generator derives the county scaffold from `exposome/conus_county_pm25_2005_2024.csv`, so the full templates match the actual CONUS county universe used elsewhere in this repo.

These templates cover the actual masked exports currently written to `output/final/`:
- `site_county_year_*.csv`
- `site_county_year_age_sex_*.csv`
- `site_year_age_sex_race_ethnicity_*.csv`
