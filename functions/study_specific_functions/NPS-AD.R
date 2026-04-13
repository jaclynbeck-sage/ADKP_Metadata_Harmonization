# Harmonize NPS-AD metadata
#
# Modifies the NPS-AD individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata files:
#   * syn55251012 (version 6, individual metadata) on Synapse
#
# NOTE: v4 of this metadata was used as the input to harmonization for GENESIS,
# and the harmonized data was re-uploaded as v5 of the same file. Then, I
# corrected some age data that was missed during that harmonization and uploaded
# it as v6, so this data is mostly harmonized already. I verified that the output
# of harmonizing v4 vs v6 are identical except:
#   1. (now) correctly-censored ages,
#   2. some geneticAncestry/geneticAncestry_isHispanic values in v6 that were
#      filled in from other data sets,
#   3. one corrected amyCerad/amyAny value (to match corrected MSBB data)
#   4. two corrected PMI values (to match corrected MSBB data)
#
# NOTE: NPS-AD determined race and ethnicity values algorithmically. To maintain
# consistency between this data set and other harmonized data sets, these values
# have been moved to new "geneticAncestry" and "geneticEthnicity" columns, and
# the original race and isHispanic columns are filled in with self-report data
# where available from other data sets.
#
# Modifications needed for version 6:
#   * Rename columns:
#     * `ethnicity` => `geneticAncestry_isHispanic`
#     * `race` => `geneticAncestry`
#   * Fix `cohort` and `dataContributionGroup` values to match data dictionary
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
