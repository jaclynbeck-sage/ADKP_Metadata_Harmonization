# Harmonize MC_snRNA
#
# Modifies the MC_snRNA individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata file: syn31563038 (version 1) on Synapse.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Change `isHispanic` values from "Caucasian" to "False"
#   * Round numerical `Braak` values down and convert to Roman numerals
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
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
harmonize_MC_snRNA <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      # Change "90_or_over" to "90+"
      ageDeath = censor_ages(ageDeath, spec),

      # There are no Hispanic individuals
      isHispanic = spec$isHispanic$hisp_false,

      # Change Braak to Roman numerals. This data set has several non-integer
      # Braak values (e.g. 0.5, 4.5), which are rounded down before conversion.
      Braak = to_Braak_stage(floor(Braak), spec),

      # Add cohort and contribution group
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo
    )
}
