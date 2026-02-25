# Harmonize MC-BrAD
#
# Modifies the MC-BrAD individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata file: syn51401700 (version 2) on Synapse.
#
# NOTE: For this function's name, "-" has been converted to "_" in the study
# name to conform to R-safe function names.
#
# NOTE: This data set has some sample overlap with AMP-AD 1.0 Mayo metadata and
# Diverse Cohorts metadata. There are some missing values in this data set that
# exist in these data sets, so those get pulled in in the main harmonization
# function.
#
# Modifications needed for version 2:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal` => `amyThal`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Round numerical `Braak` values down and convert to Roman numerals
#   * Convert `amyThal` values to conform to data dictionary
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
harmonize_MC_BrAD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal
    ) |>
    mutate(
      # Change "90_or_over" to "90+"
      ageDeath = censor_ages(ageDeath, spec),

      # Change 0 and 1 Thal values to "None" or "Phase 1"
      amyThal = case_match(amyThal,
                           0 ~ spec$amyThal$none,
                           1 ~ spec$amyThal$phase1,
                           .default = as.character(amyThal)
      ),

      # Change Braak to Roman numerals. This data set has several non-integer
      # Braak values (e.g. 0.5, 4.5), which are rounded down before conversion.
      Braak = to_Braak_stage(floor(Braak), spec),

      # Add cohort and contribution group
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo
    )
}
