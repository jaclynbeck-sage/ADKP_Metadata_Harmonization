# Harmonize NPS-AD metadata
#
# Modifies the NPS-AD individual metadata file to conform to the ADKP data
# dictionary. A separate neuropathology file needs to be merged with the
# individual metadata before harmonization to obtain Braak scores.
#
# Source metadata files:
#   * syn55251012 (version 4, individual metadata) on Synapse
#   * syn55251003 (version 1, neuropathology data) on Synapse
#
# NOTE: Cerad mapping is defined in the NPS-AD data dictionary (syn57373364).
#
# NOTE: NPS-AD determined race and ethnicity values algorithmically. To maintain
# consistency between this data set and other harmonized data sets, these values
# have been moved to new "geneticAncestry" and "geneticEthnicity" columns, and
# the original race and isHispanic columns are filled in with self-report data
# where available from other data sets.
#
# Modifications needed for version 4:
#   * Merge neuropathology data with individual metadata
#   * Rename columns:
#     * `ethnicity` => `geneticAncestry_isHispanic`
#     * `CERAD` => `amyCerad`
#     * `BRAAK_AD` => `Braak`
#     * `race` => `geneticAncestry`
#   * Remove ageDeath values of "89+" so they get replaced with the correct
#     values from Diverse Cohorts during de-duplication
#   * Fix two values in `diverseCohortsIndividualIDFormat` to match the correct
#     format in Diverse Cohorts
#   * Convert `PMI` values from minutes to hours
#   * Convert `geneticAncestry_isHispanic` values to "True" or "False"
#   * Convert Braak numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Fix `cohort` values to match data dictionary
#   * Add the `dataContributionGroup` column with values appropriate to each
#     cohort.
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
harmonize_NPS_AD <- function(metadata, neuropath, spec) {
  metadata <- metadata |>
    # Braak is all NA in the individual metadata file but has values in the
    # neuropath file, CERAD overlaps in both files so we get rid of one of them
    select(-Braak, -CERAD) |>
    merge(neuropath)

  metadata |>
    select(-Component) |>
    dplyr::rename(
      geneticAncestry = race,
      geneticAncestry_isHispanic = ethnicity,
      amyCerad = CERAD,
      Braak = BRAAK_AD
    ) |>
    mutate(
      # Remove "89+"
      ageDeath = ifelse(ageDeath == "89+", NA, ageDeath),

      # Fix to allow comparison to Diverse Cohorts
      diverseCohortsIndividualIDFormat = case_match(
        diverseCohortsIndividualIDFormat,
        29637 ~ "29637_MSSM",
        29582 ~ "29582_MSSM",
        .default = as.character(diverseCohortsIndividualIDFormat)
      ),

      # PMI from minutes to hours
      PMI = PMI / 60,
      pmiUnits = "hours",

      # Cerad mapping per NPS-AD data dictionary (syn57373364)
      amyCerad = case_match(
        amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$sparse,
        3 ~ spec$amyCerad$moderate,
        4 ~ spec$amyCerad$frequent,
        .default = as.character(amyCerad)
      ),

      # Convert Braak to Roman numerals
      Braak = to_Braak_stage(Braak, spec),

      # Change MSBB cohort values to match data dictionary. All other cohort
      # values are correct as-is, except for ROSMAP. The correct ROSMAP cohort
      # values will be filled in during de-duplication.
      cohort = ifelse(cohort == "MSBB", spec$cohort$msbb, cohort),

      # Add contribution group based on cohort values
      dataContributionGroup = case_match(
        cohort,
        spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
        spec$cohort$hbcc ~ spec$dataContributionGroup$nimh,
        c(spec$cohort$ros, spec$cohort$map) ~ spec$dataContributionGroup$rush,
        "ROSMAP" ~ spec$dataContributionGroup$rush,
        .default = ""
      )
    )
}
