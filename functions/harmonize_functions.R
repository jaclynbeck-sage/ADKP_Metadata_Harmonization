# This file contains the generic harmonize function that is called for all data
# sets. It also loads necessary libraries and sources the processing scripts
# needed to harmonize each data set. This file also contains generic value
# conversion functions for specific columns like converting Braak scores to
# Roman numerals, and functions for derived columns like amyA and bScore.

library(dplyr)
library(stringr)

# Load all data set file functions
dataset_files <- list.files(file.path("functions", "study_specific_functions"),
                            pattern = "\\.R$", full.names = TRUE)
invisible(sapply(dataset_files, source, .GlobalEnv))

# Generic harmonization function
#
# Runs the appropriate dataset-specific function to rename and harmonize
# variables, fills in any NA values with "missing or unknown", and adds any
# missing columns to the data frame.
#
# Arguments:
#   study_obj - a single study object from the "studies" list in `config.yml`
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `config.yml` file
#   extra_metadata - a `data.frame` of extra metadata that is needed by several
#     studies, which has different information depending on study. Different
#     studies have extra metadata in different file formats and it may not be
#     from Synapse, so the data frame is passed in rather than using the
#     extra_metadata field of study_obj and reading it in this function. If the
#     study does not need extra metadata, this should be NULL.
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the data dictionary.
#   Columns not defined in the data dictionary are left as-is.
#
harmonize <- function(study_obj, spec, extra_metadata = NULL) {
  # Download the main metadata file from Synapse
  meta_file <- synapse_download(study_obj$metadata_id, study_obj$name)
  metadata <- read.csv(meta_file$path)

  # Find the study-specific harmonization function to call.
  study_fn <- match.fun(study_obj$harmonize_fn)

  if (is.null(extra_metadata)) {
    metadata <- study_fn(metadata, spec)
  } else {
    metadata <- study_fn(metadata, extra_metadata, spec)
  }

  # Add any missing fields
  missing_fields <- setdiff(spec$required_columns, colnames(metadata))
  for (field in missing_fields) {
    metadata[, field] <- spec$missing
  }

  # Don't fill NA values with "missing" in the ageDeath or PMI columns
  cols_fill <- setdiff(spec$required_columns, c("ageDeath", "PMI"))

  metadata <- metadata |>
    mutate(
      # Fix fields that might be read in as numeric but should be characters
      across(any_of(cols_fill), as.character),

      # Fill NAs in character columns as "missing or unknown"
      across(any_of(cols_fill), ~ ifelse(is.na(.x), spec$missing, .x)),

      # Add or update derived columns
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyAny = get_amyAny(amyCerad, spec),
      amyA = get_amyA(amyThal, spec),
      bScore = get_bScore(Braak, spec),

      # Species should be "Human" for all studies
      species = spec$species,

      # Add study name (this variable is temporary and does not appear in the
      # final output)
      study = study_obj$name,

      # Add the original filename to keep track of it (this variable is
      # temporary and does not appear in the final output)
      filename = meta_file$name
    )

  # Put harmonized fields first in the data frame
  metadata <- metadata |>
    select(all_of(spec$required_columns), !all_of(spec$required_columns),
           # Remove this column if it exists, as it's replaced by "cohort"
           -any_of("individualIdSource"))

  return(metadata)
}


# Convert numbers to Braak stage values
#
# This function converts numerical Braak values (0-6) to "Stage " + the Roman
# numeral version of the number, as defined by the ADKP data dictionary.
#
# Arguments:
#   num - a vector containing numerical Braak stages from 0 to 6 or `NA`
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector with all Braak values converted to "Stage " + Roman numeral, and
#   all `NA` values converted to "missing or unknown". This function always
#   defaults to returning the original value as a character string if it doesn't
#   meet any of the criteria in the `case_match` statement, so the value will
#   fail validation and can be examined.
to_Braak_stage <- function(num, spec) {
  case_match(
    num,
    0 ~ spec$Braak$none,
    1 ~ spec$Braak$stage1,
    2 ~ spec$Braak$stage2,
    3 ~ spec$Braak$stage3,
    4 ~ spec$Braak$stage4,
    5 ~ spec$Braak$stage5,
    6 ~ spec$Braak$stage6,
    NA ~ spec$missing,
    .default = as.character(num)
  )
}


# Get bScore values based on Braak
#
# This function turns (harmonized) Braak scores into the corresponding values
# for bScore in the data dictionary:
#   "None" => "None"
#   "Stage I", "Stage II" => "Low (Stage I-II)"
#   "Stage III", "Stage IV" => "Moderate (Stage III-IV)"
#   "Stage V", "Stage VI" => "High (Stage V-VI)"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   Braak - a vector containing harmonized, data dictionary-compliant Braak
#           scores, which are either "Stage " + a Roman numeral or "missing or
#           unknown".
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector of strings with bScore values derived from Braak. Values should be
#   as described above.
#
#   This function always defaults to returning the original Braak value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_bScore <- function(Braak, spec) {
  case_match(
    Braak,
    spec$Braak$none ~ spec$bScore$none,
    c(spec$Braak$stage1, spec$Braak$stage2) ~ spec$bScore$stage1_2,
    c(spec$Braak$stage3, spec$Braak$stage4) ~ spec$bScore$stage3_4,
    c(spec$Braak$stage5, spec$Braak$stage6) ~ spec$bScore$stage5_6,
    NA ~ spec$missing,
    .default = as.character(Braak)
  )
}


# Get amyAny values based on amyCerad
#
# This function turns (harmonized) amyCerad scores into the corresponding values
# for amyAny in the data dictionary:
#   "None/No AD/C0" => "no"
#   "Sparse/Possible/C1", "Moderate/Probable/C2", "Frequent/Definite/C3" => "yes"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   amyCerad - a vector containing harmonized, data dictionary-compliant
#          amyCerad scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector of strings with amyAny values derived from amyCerad. Values should
#   be as described above.
#
#   This function always defaults to returning the original amyCerad value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_amyAny <- function(amyCerad, spec) {
  case_match(
    amyCerad,
    spec$amyCerad$none ~ spec$amyAny$amyNo,
    spec$amyCerad$sparse ~ spec$amyAny$amyYes,
    spec$amyCerad$moderate ~ spec$amyAny$amyYes,
    spec$amyCerad$frequent ~ spec$amyAny$amyYes,
    NA ~ spec$missing,
    .default = as.character(amyCerad)
  )
}


# Get amyA values based on amyThal
#
# This function turns (harmonized) amyThal scores into the corresponding values
# for amyA in the data dictionary:
#   "None" => "None"
#   "Phase 1", "Phase 2" => "Thal Phase 1 or 2"
#   "Phase 3" => "Thal Phase 3"
#   "Phase 4", "Phase 5" => "Thal Phase 4 or 5"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   amyThal - a vector containing harmonized, data dictionary-compliant
#          amyThal scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector of strings with amyA values derived from amyThal, with values as
#   described above.
#
#   This function always defaults to returning the original amyThal value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_amyA <- function(amyThal, spec) {
  case_match(
    amyThal,
    spec$amyThal$none ~ spec$amyA$none,
    c(spec$amyThal$phase1, spec$amyThal$phase2) ~ spec$amyA$phase1_2,
    spec$amyThal$phase3 ~ spec$amyA$phase3,
    c(spec$amyThal$phase4, spec$amyThal$phase5) ~ spec$amyA$phase4_5,
    NA ~ spec$missing,
    .default = as.character(amyThal)
  )
}


# Get APOE4 status based on genotype
#
# This function turns (harmonized) apoeGenotype scores into the corresponding
# values for apoe4Status in the data dictionary:
#   "22", "23", "33" => "no"
#   "24", "34", "44" => "yes"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   apoeGenotype - a vector containing harmonized, data dictionary-compliant
#          apoeGenotype scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector of strings with apoe4Status values derived from apoeGenotype, with
#   values as described above.
#
#   This function always defaults to returning the original apoeGenotype value
#   as a character string if it doesn't meet any of the criteria in the
#   `case_match` statement, so the value will fail validation and can be
#   examined.
get_apoe4Status <- function(apoeGenotype, spec) {
  case_match(
    apoeGenotype,
    # "yes" if there is a 4
    c(
      spec$apoeGenotype$e2e4,
      spec$apoeGenotype$e3e4,
      spec$apoeGenotype$e4e4
    ) ~ spec$apoe4Status$e4yes,
    # "no" if there is no 4
    c(
      spec$apoeGenotype$e2e2,
      spec$apoeGenotype$e2e3,
      spec$apoeGenotype$e3e3
    ) ~ spec$apoe4Status$e4no,
    NA ~ spec$missing,
    .default = as.character(apoeGenotype)
  )
}


# Censor ages 90 or above
#
# This function censors any ages 90 or above by replacing them with "90+". It
# also converts some studies' versions of this from "90_or_over" to "90+". Empty
# strings and "missing or unknown" values should be replaced with `NA`.
#
# Arguments:
#   ages - a vector of ages, which may be strings or numerical
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a vector of strings with age values properly censored
censor_ages <- function(ages, spec) {
  case_when(
    ages %in% c("90+", "90_or_over") ~ spec$ageDeath$over90,
    ages == "" ~ NA,
    ages == spec$missing ~ NA,
    suppressWarnings(as.numeric(ages)) >= 90 ~ spec$ageDeath$over90,
    .default = as.character(ages)
  )
}


# Determine the value of ADoutcome for Diverse Cohorts data
#
# ADoutcome is not changed for data not from Diverse Cohorts, and data from
# Diverse Cohorts but contributed by Mayo.
#
# Arguments:
#   .data - a single row of a data frame, as from inside a rowwise() |> mutate() statement
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   a single value for ADoutcome, one of "Control", "AD", "Other", or "missing or unknown"
determineADoutcome <- function(.data, spec) {
  # Don't make an ADoutcome value for non-Diverse Cohorts data, and don't change
  # ADoutcome for samples coming from Mayo Clinic
  if (.data[["study"]] != "AMP-AD_DiverseCohorts" |
      .data[["derivedOutcomeBasedOnMayoDx"]] == TRUE) {
    return(.data[["ADoutcome"]])
  }

  high_Braak <- .data[["Braak"]] %in% c(spec$Braak$stage4,
                                        spec$Braak$stage5,
                                        spec$Braak$stage6)
  low_Braak <- .data[["Braak"]] %in% c(spec$Braak$none,
                                       spec$Braak$stage1,
                                       spec$Braak$stage2,
                                       spec$Braak$stage3)

  high_amyCerad <- .data[["amyCerad"]] %in% c(spec$amyCerad$moderate,
                                              spec$amyCerad$frequent)
  low_amyCerad <- .data[["amyCerad"]] %in% c(spec$amyCerad$none,
                                             spec$amyCerad$sparse)

  ADoutcome <- case_when(
    # AD: Braak IV-VI & Cerad Moderate or Frequent
    high_Braak & high_amyCerad ~ "AD",

    # Control: Braak 0-III & Cerad None or Sparse
    low_Braak & low_amyCerad ~ "Control",

    # Other: Braak IV-VI & Cerad None or Sparse
    high_Braak & low_amyCerad ~ "Other",

    # Other: Braak 0-III & Cerad Moderate or Frequent
    low_Braak & high_amyCerad ~ "Other",

    # Missing one or both of Braak and Cerad
    .data[["Braak"]] == spec$missing |
      .data[["amyCerad"]] == spec$missing ~ spec$missing,

    .default = .data[["ADoutcome"]]
  )

  return(ADoutcome)
}


# Update ADoutcome values for Diverse Cohorts data in case Braak or amyCerad
# values changed during de-duplication.
#
# ADoutcome is not changed for data not from Diverse Cohorts, and data from
# Diverse Cohorts but contributed by Mayo.
#
# Arguments:
#   meta_all - dataframe of harmonized data, which may contain data from any study
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#
# Returns:
#   meta_all with ADoutcome values updated for Diverse Cohorts data only
updateADoutcome <- function(meta_all, spec) {
  meta_all |>
    rowwise() |>
    mutate(ADoutcome = determineADoutcome(.data, spec)) |>
    ungroup() |>
    as.data.frame()
}
