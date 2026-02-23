# Harmonize MayoRNAseq metadata
#
# Modifies the original MayoRNAseq individual metadata to conform to the ADKP
# data dictionary.
#
# Source metadata file: syn23277389 (version 7) on Synapse.
#
# Modifications needed for version 7:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal` => `amyThal`
#   * Change `ageDeath` value "90_or_over" to "90+"
#   * Rename `isHispanic` value "Caucasian" to "False"
#   * Replace `amyThal` value "0" with "None", and add "Phase" in front of other
#       `amyThal` numerical values
#   * Convert `Braak` value "0" with "None", and convert other numerical `Braak`
#       values to "Stage " + a Roman numeral
#   * Add columns `dataContributionGroup` = "Mayo" and `cohort` = "Mayo Clinic"
#
# NOTE: There is one individual with an incorrect `race` value, which is
# corrected manually here based on updated information from Diverse Cohorts.
#
# NOTE: all amyCerad values are NA in this data set and will get filled in with
# "missing or unknown" in the main harmonize() function.
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
harmonize_MayoRNAseq <- function(metadata, biospecimen, spec) {
  cohort_info <- metadata |>
    select(individualID) |>
    merge(biospecimen) |>
    select(individualID, specimenIdSource) |>
    distinct() |>
    subset(!is.na(specimenIdSource) & specimenIdSource != "SNPs") |>
    dplyr::rename(cohort = specimenIdSource)

  metadata |>
    merge(cohort_info, all.x = TRUE, sort = FALSE) |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal
    ) |>
    mutate(
      # Change "90_or_over" to "90+"
      ageDeath = censor_ages(ageDeath, spec),

      # Change "Caucasian" to False
      isHispanic = ifelse(isHispanic == "Caucasian",
                          spec$isHispanic$hisp_false,
                          isHispanic
      ),

      ## Manual correction
      race = ifelse(individualID == "11387", spec$race$other, race),
      ##

      # Add "Phase" in front of Thal values
      amyThal = case_match(amyThal,
                           NA ~ spec$missing,
                           0 ~ spec$amyThal$none,
                           .default = paste("Phase", amyThal)
      ),

      # Change Braak to Roman numerals. This data set has several non-integer
      # Braak values (e.g. 0.5, 4.5), which are rounded down before conversion.
      Braak = to_Braak_stage(floor(Braak), spec),

      # Map cohort values from biospecimen file
      cohort = case_match(
        cohort,
        "BannerSun" ~ spec$cohort$banner,
        "MayoBrainBank" ~ spec$cohort$mayo,
        "UniversityKentucky" ~ spec$cohort$ukentucky,
        NA ~ spec$cohort$mayo,
        .default = cohort
      ),

      # Add contribution group
      dataContributionGroup = spec$dataContributionGroup$mayo,
    )
}
