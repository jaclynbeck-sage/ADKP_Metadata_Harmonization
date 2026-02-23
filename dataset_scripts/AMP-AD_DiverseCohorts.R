# Harmonize Diverse Cohorts metadata
#
# Makes minor edits to the Diverse Cohorts individual metadata file, which was
# previously harmonized using a very similar data dictionary to what is needed
# for the ADKP harmonization project.
#
# Source metadata file: syn51757646 (version 21) on Synapse.
#
# NOTE: For this function's name, "-" has been converted to "_" in the study
# name to conform to R-safe function names.
#
# Modifications needed for version 21:
#   * Change `ageDeath` and `PMI` value "missing or unknown" to `NA`
#   * Change `PMI` to a numeric column
#   * Change `race` "other" to "Other"
#   * Rename `isHispanic` values from ["TRUE", "FALSE"] to ["True", "False"]
#   * Update `dataContributionGroup` and `cohort` to conform to the data dictionary
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `config.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the data dictionary.
#   Columns not defined in the data dictionary are left as-is.
#
harmonize_AMP_AD_DiverseCohorts <- function(metadata, spec) {
  metadata |>
    dplyr::mutate(
      # Ages should already be censored, but this also changes "missing" to NA
      ageDeath = censor_ages(ageDeath, spec),

      # Change PMI to numeric
      PMI = ifelse(PMI == spec$missing, NA,
                   suppressWarnings(as.numeric(PMI))
      ),

      # Change "other" to "Other"
      race = ifelse(race == "other", "Other", race),

      # Change isHispanic from all caps to title case
      isHispanic = case_match(
        as.character(isHispanic),
        "TRUE" ~ spec$isHispanic$hisp_true,
        "FALSE" ~ spec$isHispanic$hisp_false,
        .default = as.character(isHispanic)
      ),

      # Standardize data contribution group
      dataContributionGroup = case_match(
        dataContributionGroup,
        "Columbia" ~ spec$dataContributionGroup$columbia,
        "MSSM" ~ spec$dataContributionGroup$mssm,
        "Rush" ~ spec$dataContributionGroup$rush,
        "Emory" ~ spec$dataContributionGroup$emory,
        "Mayo" ~ spec$dataContributionGroup$mayo,
        .default = dataContributionGroup
      ),

      # Standardize cohort values
      cohort = case_match(
        cohort,
        "Banner" ~ spec$cohort$banner,
        "Mayo Clinic" ~ spec$cohort$mayo,
        "Mt Sinai Brain Bank" ~ spec$cohort$msbb,
        .default = cohort
      )
    )
}
