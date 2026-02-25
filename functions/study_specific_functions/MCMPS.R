# Harmonize MCMPS
#
# Modifies the MCMPS individual metadata file to conform to the ADKP data
# dictionary.
#
# Source metadata file: syn25891193 (version 1) on Synapse.
#
# NOTE: These samples come from living tissue and do not have a value for
# `ageDeath` or `PMI`.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Trim extra spaces from `race` values
#   * Trim extra spaces from `isHispanic` values and convert to True/False
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
harmonize_MCMPS <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      # Remove extra spaces in race values
      race = str_trim(race),

      # Remove extra spaces in isHispanic values and convert to True/False
      isHispanic = str_trim(isHispanic),
      isHispanic = case_match(isHispanic,
                              "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
                              "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
                              "Middle Eastern" ~ spec$isHispanic$hisp_false,
                              .default = isHispanic
      ),

      # Add cohort and dataContributionGroup
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo
    )
}
