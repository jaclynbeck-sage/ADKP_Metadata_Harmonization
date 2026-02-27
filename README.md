# ADKP Metadata Harmonization
This code harmonizes metadata from studies across the [AD Knowledge Portal](https://adknowledgeportal.synapse.org/). Common columns between studies (for example age, sex, Braak stage) have been standardized to have the same column names and values across all studies to make cross-comparison simpler. The output of this project is the [ADKP_Harmonization_Study](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn70781457). A full study description, a list of included studies, and a methods overview are available on the study page. 

The code has very detailed documentation in every file, so this README contains a general overview of each step plus where to look for specific implementation details.

## Code overview

The general process is as follows:
1. Harmonize each study
2. Find overlapping individuals between studies
3. Compare data for overlapping individuals and fill in missing information where available
4. Write the harmonized data to files and upload them to the project on Synapse

Due to the scope of this project, the process is broken out into several sub-functions across multiple files:
* `config.yml` - a configuration file that defines study names, study file IDs, and column requirements
* `docs/ADKP_harmonized_data_dictionary.docx` - an editable version of the data dictionary
* `Harmonize_All_Studies.R` - the top-level or "main" file
* `functions` - all sub-functions are in this folder
  * `deduplication_functions.R` - functions to handle filling in missing values between overlapping individuals
  * `harmonize_functions.R` - generic harmonization operations that are not study-specific
  * `overlap_functions.R` - functions for finding overlapping ID between studies and verifying they refer to the same individual
  * `util_functions.R` - functions for file I/O, Synapse I/O, and validation of the final harmonized files
  * `study_specific_functions/*` - one file/function for each study that handles study-specific changes needed for harmonization

----
### 1. Harmonize each study

There is a generic `harmonize` function (in `functions/harmonize_functions.R`) that downloads the main metadata file for a given study, calls that study's specific harmonization function (located in `functions/study_specific_functions`), and then does the following final changes:
1. Add any missing columns and fill them with "missing or unknown"
2. Change any `NA` values in character columns to "missing or unknown"
3. Fill in values for "derived" columns. These are columns like `amyA` or `bScore`, which derive their value from another column (`amyThal`, `Braak`) rather than being provided by the data contributor.
4. Put all the "harmonized" columns first in the data frame, and put the non-harmonized columns at the end

Each study-specific harmonization file covers any of the following as needed for a given study:
1. Column renames
2. Value re-mapping
3. Spelling, capitalization, or punctuation changes
4. Censoring ages above 90
5. Creation or modification of `dataContributionGroup` and `cohort` values

----
### 2. Find overlapping individuals

Many of the studies we harmonized share some of the same individuals, however their `individualIDs` may be slightly differently-formatted and/or some study data is more complete than others. We do *not* change `individualIDs`, however we do try to fill in missing information if at least one study has non-missing data in that field for a shared individual. The process to find which studies share each individual is as follows:
1. *Temporarily* alter all `individualIDs` to be in the same format for direct comparison. The only major alterations necessary are for MSBB IDs in MSBB, NPS-AD, and Diverse Cohorts data, which have the same core numerical ID but in different string formats.
2. Find all pairs of studies that have the same ID
3. For every pair, get the individual's row of values from each study's data frame and compare the data in overlapping columns.
    * The data has been harmonized by this point, so this should include all harmonized columns plus any non-harmonized columns that happen to have the same column name.
    * ID columns and derived columns (which are redundant with their source columns) are excluded from comparison.
4. Count up the number of matching values, the number of columns where one or both studies is missing a value, and the number of non-missing mismatched values.
5. Determine which studies with the same ID actually represent the same individual. Studies share the same individual if they have no mis-matched values and less than 50% of their columns are missing one or both values.
    * Some special case individuals from NPS-AD and ROSMAP that do not meet these criteria were verified as matching by hand and added to the list of true overlaps.
    * In a small number of cases, two studies had the same ID but for clearly different individuals. These cases were all verified by hand as non-matching.
6. Create a "group ID" representing each individual, and create a data frame that maps it to the `individualID` from each study that contains that individual. Write this to a file.
    * Some individuals are shared by more than two studies, so each group ID maps to *all* studies that share the same individual, not just pairs of studies.
  
This process is entirely contained in the `functions/overlap_functions.R` file.

----
### 3. Fill in missing information

For each individual (represented by their group ID), we compared the data across all studies that contain that individual. If one study had missing data in a field and at least one other study had a non-missing value, we filled in the missing value with the non-missing one. If two or more studies shared a non-missing value, we confirmed that both values matched before replacing missing values. 

**NOTE:** There are a few special cases (3 individuals total) where there are still mis-matching values between NPS-AD data and overlapping studies. These have been left as-is until the corrected values can be obtained from the data contributor.  

This process takes place in the `functions/deduplication_functions.R` file. 

----
### 4. Write to files and upload

The final harmonized files contain all harmonized columns *plus* any non-harmonized columns that were in the original metadata for that study. We do not add any non-harmonized columns from other data sets, even when there is a large number of shared individuals. 

A harmonized file for each study, plus the data dictionary and the data frame mapping shared individuals across studies, are available [on the ADKP/Synapse](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn70781457).
