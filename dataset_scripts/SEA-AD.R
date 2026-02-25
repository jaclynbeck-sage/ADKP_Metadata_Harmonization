# Harmonize SEA-AD metadata
#
# Modifies the SEA-AD individual metadata file to conform to the ADKP data
# dictionary. The version of SEA-AD that is on Synapse is missing
# Hispanic/Latino information that is present in the version released by the
# Allen Institute on brain-map.org. We use the version on Synapse, because it
# has been curated and approved for release in the AD Knowledge Portal, but pull
# in the missing Hispanic/Latino information from the Allen Institute version.
#
# Source metadata files:
#   * syn31149116 (version 8) on Synapse
#   * the "Donor Metadata" file downloaded from https://portal.brain-map.org.
#     Full URL: https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx
#
# NOTE: For this function's name, "-" has been converted to "_" in the study
# name to conform to R-safe function names.
#
# Modifications needed for version 8:
#   * Rename several columns:
#     * `pmi` => `PMI`
#     * `Hispanic.Latino` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal phase` => `amyThal`
#   * Convert `isHispanic` values ["Yes", "No"] to ["True", "False"]
#   * Individuals may have multiple races in the `race` column, separated by a
#       semi-colon. Any value with "American Indian" in it is assigned to
#       "American Indian or Alaska Native".
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyThal` values from "Thal #" to "None" or "Phase #"
#   * Add `cohort` ("SEA-AD"), `dataContributionGroup` ("Allen Institute")
#       and `study` ("SEA-AD")
#   * Use Atherosclerosis values from the Allen metadata (even though these
#     aren't part of the harmonization), because v8 of the Synapse metadata
#     incorrectly combines "None" and NA values.
#
# Arguments:
#   metadata_synapse - a `data.frame` of metadata from the Synapse metadata
#     file. Columns are variables and rows are individuals.
#   metadata_allen - a `data.frame` of metadata from brain-map.org. Columns are
#     variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `config.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the data dictionary.
#   Columns not defined in the data dictionary are left as-is.
#
harmonize_SEA_AD <- function(metadata_synapse, metadata_allen, spec) {
  # Keep only the IDs, Hispanic/Latino column, and Atherosclerosis column from
  # the Allen Institute data
  metadata_allen <- metadata_allen |>
    select("Donor ID", "Hispanic/Latino", "Atherosclerosis")

  metadata_synapse <- metadata_synapse |> select(-Atherosclerosis, -X) # v8 fix

  # v8 fix: match column names of v7
  # TODO is this necessary?
  colnames(metadata_synapse) <- colnames(metadata_synapse) |>
    str_replace_all("\\.", " ")

  meta <- merge(metadata_synapse, metadata_allen,
                by.x = "individualID",
                by.y = "Donor ID",
                all = TRUE
  ) |>
    dplyr::relocate(Atherosclerosis, .before = Arteriolosclerosis) # v8 fix

  # v8 fix: Set all empty strings to NA to match previous versions
  meta[meta == ""] <- NA

  meta |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = "Hispanic/Latino",
      amyCerad = CERAD,
      amyThal = "Thal phase",
      `LATE-NC stage` = "LATE NC stage" # v8 fix to match v7
    ) |>
    dplyr::mutate(
      # Convert Yes/No to True/False
      isHispanic = case_match(isHispanic,
                              "No" ~ spec$isHispanic$hisp_false,
                              "Yes" ~ spec$isHispanic$hisp_true,
                              "Unknown" ~ spec$missing,
                              .default = isHispanic
      ),

      # Assign any multi-race values with "American Indian" to the American
      # Indian category
      race = case_when(
        grepl("American Indian", race) ~ spec$race$Amer_Ind,
        .default = race
      ),

      # Re-map amyCerad values
      amyCerad = case_match(amyCerad,
                            0 ~ spec$amyCerad$none,
                            1 ~ spec$amyCerad$sparse,
                            2 ~ spec$amyCerad$moderate,
                            3 ~ spec$amyCerad$frequent,
                            .default = as.character(amyCerad)
      ),

      # Change "Thal X" to "Phase X"
      amyThal = ifelse(amyThal == "Thal 0",
                       spec$amyThal$none,
                       str_replace(amyThal, "Thal", "Phase")
      ),

      # Convert Braak to Roman numerals
      Braak = to_Braak_stage(Braak, spec),

      # Add cohort and contribution group
      cohort = spec$cohort$sea_ad,
      dataContributionGroup = spec$dataContributionGroup$allen_institute,
    )
}
