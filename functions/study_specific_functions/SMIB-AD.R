# Harmonize SMIB-AD
#
# Modifies the SMIB-AD individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata file: syn22432749 (version 1) on Synapse.
#
# NOTE: CERAD values are coded using ROSMAP's system.
#
# NOTE: For this function's name, "-" has been converted to "_" in the study
# name to conform to R-safe function names.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Censor `ageDeath` values over 90
#   * Change `race` values from "European" to "White"
#   * Change `isHispanic` values from "European" to "False"
#   * Convert Braak numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add a `cohort` column with either "SMRI" or "Banner"
#   * Add a `dataContributionGroup` column with values for Stanley and Banner
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
harmonize_SMIB_AD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      # Censor ages over 90
      ageDeath = censor_ages(ageDeath, spec),

      # Change "European" to "White"
      race = spec$race$White,

      # There are no Hispanic individuals
      isHispanic = spec$isHispanic$hisp_false,

      # Convert Braak to Roman numerals
      Braak = to_Braak_stage(Braak, spec),

      # Re-map amyCerad
      amyCerad = case_match(amyCerad,
                            1 ~ spec$amyCerad$none,
                            2 ~ spec$amyCerad$sparse,
                            4 ~ spec$amyCerad$frequent,
                            .default = as.character(amyCerad)
      ),

      # Add cohort and contribution group
      cohort = ifelse(is.na(individualIdSource),
                      spec$cohort$smri,
                      spec$cohort$banner
      ),
      dataContributionGroup = ifelse(cohort == spec$cohort$smri,
                                     spec$dataContributionGroup$stanley,
                                     spec$dataContributionGroup$banner
      )
    )
}
