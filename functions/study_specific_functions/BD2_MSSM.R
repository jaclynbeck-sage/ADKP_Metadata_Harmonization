# Harmonize BD2_MSSM metadata
#
# Makes minor edits to the BD2_MSSM individual metadata file.
#
# Source metadata file: syn76211469 (version 19) on Synapse.
#
# Modifications needed for version 19:
#   * Rename multiple columns:
#     * `ethnicity` => `isHispanic`
#     * `bScore` => `bScore_NFT`
#   * Change `race` "other" to "Other"
#   * Update `dataContributionGroup` and `cohort` to conform to the data dictionary
#   * Get `isHispanic` values from `ethnicity`
#   * Change `ageDeath` and `PMI` value "missing or unknown" to `NA`
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Remove redundant `Braak_NFT` column
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
harmonize_BD2_MSSM <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      bScore = bScore_NFT,
      isHispanic = ethnicity
    ) |>

    select(-Braak_NFT) |>

    dplyr::mutate(
      # Change "other" to "Other"
      race = ifelse(race == "other", spec$race$other, race),

      # Ages should already be censored, but this also changes "missing" to NA
      ageDeath = censor_ages(ageDeath, spec),

      # Change PMI to numeric
      PMI = suppressWarnings(as.numeric(PMI)),

      # Get isHispanic from ethnicity
      isHispanic = case_match(
        isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
        "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
        "unknown" ~ spec$missing,
        .default = isHispanic
      ),

      # Make Braak stages into Roman numerals. Unharmonized Braak values are
      # in the format "BRAAK Stage 1", "No BRAAK stage", etc.
      Braak = to_Braak_stage(
        suppressWarnings(
          str_replace(Braak, ".* ", "") |> as.numeric()
        ), spec
      ),

      # Standardize data contribution group
      dataContributionGroup = case_match(
        dataContributionGroup,
        "Icahn School of Medicine at Mount Sinai" ~ spec$dataContributionGroup$mssm,
        .default = dataContributionGroup
      ),

      # Standardize cohort
      cohort = ifelse(cohort == "UPitt", spec$cohort$upitt_nbb, cohort)
    )
}
