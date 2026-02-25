# Harmonize ROSMAP metadata
#
# Modifies the ROSMAP individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata file: syn3191087 (version 11) on Synapse.
#
# The clinical codebook describing the meanings of values in each column can be
# found at syn3191090 on Synapse.
#
# Modifications needed for version 11:
#   * Rename multiple columns:
#     * `spanish` => `isHispanic`
#     * `age_death` => `ageDeath`
#     * `pmi` => `PMI`
#     * `msex` => `sex`
#     * `apoe_genotype` => `apoeGenotype`
#     * `ceradsc` => `amyCerad`
#     * `braaksc` => `Braak`
#     * `Study` => `cohort`
#   * Convert `ageDeath` empty string values to `NA`
#   * Convert `sex` values [0, 1] to ["female", "male"]
#   * Convert `race` numerical values to values in data dictionary
#   * Convert `isHispanic` values [1, 2] to ["True", "False"]
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add `dataContributionGroup` = "Rush" and `study` = "ROSMAP"
#
# NOTE: There is one individual with an incorrect value for `isHispanic`, which
# is manually corrected here based on updated data from Diverse Cohorts.
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
harmonize_ROSMAP <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      isHispanic = spanish,
      ageDeath = age_death,
      PMI = pmi,
      sex = msex,
      apoeGenotype = apoe_genotype,
      amyCerad = ceradsc,
      Braak = braaksc,
      cohort = Study
    ) |>

    # Change sex 0 and 1 values to "male" and "female"
    dplyr::mutate(
      sex = case_match(sex,
                       1 ~ spec$sex$male,
                       0 ~ spec$sex$female,
                       .default = as.character(sex)
      ),

      # Ages are already censored, but this also changes "" values to NA
      ageDeath = censor_ages(ageDeath, spec),

      # Translate numbers to full race names
      race = case_match(race,
                        1 ~ spec$race$White,
                        2 ~ spec$race$Black,
                        3 ~ spec$race$Amer_Ind,
                        4 ~ spec$race$Pacif_Islander,
                        5 ~ spec$race$Asian,
                        6 ~ spec$race$other,
                        7 ~ spec$missing,
                        .default = as.character(race)
      ),

      # Translate numbers to True/False
      isHispanic = case_match(isHispanic,
                              1 ~ spec$isHispanic$hisp_true,
                              2 ~ spec$isHispanic$hisp_false,
                              .default = as.character(isHispanic)
      ),
      ## manual correction
      isHispanic = ifelse(individualID == "R8412417",
                          spec$isHispanic$hisp_false,
                          isHispanic
      ),
      ##

      # Convert Braak to Roman numerals
      Braak = to_Braak_stage(Braak, spec),

      # Re-map amyCerad values
      amyCerad = case_match(amyCerad,
                            1 ~ spec$amyCerad$frequent,
                            2 ~ spec$amyCerad$moderate,
                            3 ~ spec$amyCerad$sparse,
                            4 ~ spec$amyCerad$none,
                            .default = as.character(amyCerad)
      ),

      # Add contribution group
      dataContributionGroup = spec$dataContributionGroup$rush,
    )
}
