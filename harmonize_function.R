# This file contains the generic harmonize function that is called for all data
# sets. It also loads necessary libraries and sources the processing scripts
# needed to harmonize each data set.

library(dplyr)
library(stringr)
source("util_functions.R")

# Load all data set file functions
dataset_files <- list.files("dataset_scripts", pattern = "\\.R$", full.names = TRUE)
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
