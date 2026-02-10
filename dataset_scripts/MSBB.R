# Harmonize MSBB metadata
#
# Modifies the original MSBB individual metadata to conform to the ADKP
# data dictionary.
#
# Source metadata file: syn6101474 (version 10) on Synapse.
#
# Modifications needed for version 10:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Convert `PMI` from minutes to hours
#   * Change `isHispanic` and `race` values to conform to the data dictionary
#   * Add `dataContributionGroup` = "Mount Sinai School of Medicine",
#     `cohort` = "Mount Sinai Brain Bank", and `study` = "MSBB"
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
harmonize_MSBB <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      # PMI is in minutes
      PMI = PMI / 60,

      # Change from race abbreviations to True/False
      isHispanic = case_match(isHispanic,
                              c("A", "B", "O", "W") ~ spec$isHispanic$hisp_false,
                              "H" ~ spec$isHispanic$hisp_true,
                              "U" ~ spec$missing,
                              .default = isHispanic
      ),

      # Change from race abbreviations to full names
      race = case_match(race,
                        "A" ~ spec$race$Asian,
                        "B" ~ spec$race$Black,
                        "H" ~ spec$race$other,
                        "W" ~ spec$race$White,
                        "O" ~ spec$race$other,
                        "U" ~ spec$missing,
                        .default = race
      ),

      # Re-map amyCerad values
      amyCerad = case_match(amyCerad,
                            1 ~ spec$amyCerad$none,
                            2 ~ spec$amyCerad$frequent,
                            3 ~ spec$amyCerad$moderate,
                            4 ~ spec$amyCerad$sparse,
                            .default = as.character(amyCerad)
      ),

      # Change Braak to Roman numerals
      # TODO do we need to floor for this data set?
      Braak = to_Braak_stage(floor(Braak), spec),

      # Add cohort and contribution group
      dataContributionGroup = spec$dataContributionGroup$mssm,
      cohort = spec$cohort$msbb
    )
}
