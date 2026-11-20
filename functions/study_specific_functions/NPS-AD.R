# Harmonize NPS-AD metadata
#
# Modifies the NPS-AD individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata files:
#   * syn55251012 (version 9, individual metadata) on Synapse
#
# NOTE: NPS-AD determined race and ethnicity values algorithmically. To maintain
# consistency between this data set and other harmonized data sets, these values
# have been moved to new "geneticAncestry" and "geneticEthnicity" columns, and
# the original race and isHispanic columns are filled in with self-report data
# where available from other data sets.
#
# Modifications needed for version 9:
#   * Rename columns:
#     * `ethnicity` => `geneticAncestry_isHispanic`
#     * `race` => `geneticAncestry`
#   * Fix `cohort` and `dataContributionGroup` values to match data dictionary
#   * Make `PMI` numeric (and convert "missing or unknown" to NA)
#   * Standardize Braak values
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   neuropath - a `data.frame` of neuropathology data for each individual, which
#     can be matched to `metadata` by `individualID`. Columns are variables and
#     rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `config.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the data dictionary.
#   Columns not defined in the data dictionary are left as-is.
#
harmonize_NPS_AD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      geneticAncestry = race,
      geneticAncestry_isHispanic = isHispanic
    ) |>
    mutate(
      # Ages should already be censored, but just in case.
      ageDeath = censor_ages(ageDeath, spec),

      # PMI is a string because there's a "missing or unknown" value
      PMI = suppressWarnings(as.numeric(PMI)),

      # Most Braak stages are correct but a few are "Stage #" for newly-added
      # samples in v9
      Braak = ifelse(
        grepl("[0-9]", Braak),
        to_Braak_stage(
          suppressWarnings(
            str_replace(Braak, "Stage ", "") |> as.numeric()
          ), spec),
        Braak
      ),

      # Change MSBB cohort values to match data dictionary. All other cohort
      # values are correct as-is.
      cohort = ifelse(cohort == "Mt Sinai Brain Bank",
                      spec$cohort$msbb, cohort),

      # Fix a few dataContributionGroup values
      dataContributionGroup = case_match(
        dataContributionGroup,
        "MSSM" ~ spec$dataContributionGroup$mssm,
        "Rush" ~ spec$dataContributionGroup$rush,
        .default = dataContributionGroup
      )
    )
}
