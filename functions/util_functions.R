# This file contains utility functions and variables that are used in multiple
# steps/scripts and are not dataset-specific. Functions include things like
# Synapse upload/download, quality control printouts, and harmonized value
# validation.
library(dplyr)
library(stringr)
library(synapser)

# Quality control printouts
#
# Prints summaries of all expected columns for visual inspection:
#   1. Reports which columns are missing from the data frame
#   2. Reports how many NA values are in each column
#   3. Prints tables of unique values in each column vs the count of each value,
#      excluding `ageDeath` and `PMI`.
#   4. Checks if there are any `ageDeath` values of 90 or over that are not
#      censored as "90+".
#
# Arguments:
#   df - a data.frame where rows are individuals and columns are variables
#   <X>_col - the real name of the column in `df` that corresponds to <X>.
#             Default values are provided, but many data sets have column names
#             that are different from the desired harmonized name (i.e. "CERAD"
#             instead of "amyCerad" or "ethnicity" instead of "isHispanic"). If
#             no column corresponding to <X> exists, `<X>_col` can be left as
#             the default value and will be reported as missing by this
#             function.
# Returns:
#   nothing
print_qc <- function(df, spec, column_renames = NULL) {
  required_columns <- spec$required_columns
  names(required_columns) <- required_columns

  # Sub in the column names that actually exist in the data frame for the
  # required columns
  if (!is.null(column_renames) && length(column_renames) > 0) {
    required_columns[names(column_renames)] <- as.character(column_renames)
  }

  missing <- setdiff(required_columns, colnames(df))
  cat("Missing columns:", paste(missing, collapse = ", "), "\n\n")

  cat("NA value check:\n\n")
  for (col_name in required_columns) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "has", sum(is.na(df[, col_name])), "NA values\n")
    }
  }

  cat("\nUnique value check:\n\n")
  for (col_name in required_columns) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "values:\n")
      tmp <- df |>
        group_by_at(col_name) |>
        count() |>
        ungroup() |>
        mutate_if(is.character, ~ paste0("\"", .x, "\"")) |>
        data.frame()

      # Limit printing if column has a lot of unique values, like an ID column
      # or an age column
      if (nrow(tmp) > 20) {
        cat("\tNot displayed: More than 20 unique values.")
      } else {
        print(tmp)
      }

      cat("\n\n")
    }
  }

  ageDeath_col <- required_columns["ageDeath"]
  if (ageDeath_col %in% colnames(df)) {
    cat("Age check:\n")
    ages_remove <- c("90+", "90_or_over", "Missing or unknown", "missing or unknown")
    tmp <- subset(df, !(df[, ageDeath_col] %in% ages_remove)) |>
      mutate(ageDeath = suppressWarnings(as.numeric(.data[[ageDeath_col]]))) |>
      subset(ageDeath >= 90)

    cat(ageDeath_col, "has", nrow(tmp), "uncensored ages.\n")
  }
}


# Validation of harmonized results
#
# Given a data frame of harmonized metadata, this function validates the following:
#   1. There are no values in any harmonized field that aren't in the data dictionary
#   2. There are no `ageDeath` values above 89
#   3. The `ageDeath` and `PMI` columns only have numbers, NAs, or "90+"
#   4. Columns whose values are derived from other columns (`apoe4Status`,
#      `amyA`, `amyAny`, and `bScore`) have the correctly-derived values. This
#      check is needed to catch any accidental differences introduced by filling
#      in missing data from another data set.
#
# If a column passes validation, the phrase "OK <column>" will print out if
# `verbose = TRUE`. If a column fails validation, the function will print out
# "X" plus a message describing the failure and failing values.
#
# Arguments:
#   metadata - a data frame of harmonized metadata, where rows are individuals
#          and columns are variables
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `config.yml` file
#   verbose - if TRUE, the phrase "OK <column>" will print out if the column
#          passes validation. If FALSE, nothing will print out for columns that
#          pass. Columns that fail validation always have a print out, so
#          setting verbose = FALSE will result in only failures being printed.
#
# Returns:
#   nothing
validate_values <- function(metadata, spec, verbose = TRUE) {
  # ageDeath should have only NA, numbers, or "90+". No numbers should be above
  # 89.
  ageDeath <- na.omit(metadata$ageDeath) |>
    setdiff(spec$ageDeath$over90) |>
    as.numeric()

  if (any(ageDeath >= 90, na.rm = TRUE)) {
    cat("X  ageDeath has uncensored ages above 90.\n")
  } else if (any(is.na(ageDeath))) {
    # NAs introduced by coercion to numeric. Before coercion, we removed all NAs
    # that existed in the column, so if NAs exist *after* coercion, that means
    # there was a string value that couldn't be converted to a numeric value.
    nas <- which(is.na(ageDeath))
    tmp <- na.omit(metadata$ageDeath) |>
      setdiff(spec$ageDeath$over90)
    cat(
      "X  ageDeath has invalid age values:",
      paste(tmp[nas], collapse = ", "), "\n"
    )
  } else if (verbose) {
    cat("OK ageDeath\n")
  }

  # PMI should have only NA or numbers
  pmi <- na.omit(metadata$PMI)
  if (length(pmi) > 0 && !is.numeric(pmi)) {
    cat("X  PMI is not numeric\n")
  } else if (verbose) {
    cat("OK PMI\n")
  }

  # Other columns should have only string values that exist in the dictionary
  cols_check <- setdiff(spec$required_columns, c("individualID", "ageDeath", "PMI"))

  for (col_name in cols_check) {
    values <- metadata[, col_name]
    expected <- c(spec$missing, unlist(spec[[col_name]]))
    if (!all(values %in% expected)) {
      cat(
        "X ", col_name, "has unexpected values:",
        paste(setdiff(values, expected), collapse = ", "), "\n"
      )
    } else if (verbose) {
      cat("OK", col_name, "\n")
    }
  }

  # Also check agreement between related columns
  if (!identical(metadata$apoe4Status, get_apoe4Status(metadata$apoeGenotype, spec))) {
    cat("X apoe4Status does not match apoeGenotype\n")
  } else if (verbose) {
    cat("OK apoeGenotype vs apoe4Status\n")
  }

  if (!identical(metadata$amyA, get_amyA(metadata$amyThal, spec))) {
    cat("X amyA does not match amyThal\n")
  } else if (verbose) {
    cat("OK amyThal vs amyA\n")
  }

  if (!identical(metadata$amyAny, get_amyAny(metadata$amyCerad, spec))) {
    cat("X amyAny does not match amyCerad\n")
  } else if (verbose) {
    cat("OK amyCerad vs amyAny\n")
  }

  if (!identical(metadata$bScore, get_bScore(metadata$Braak, spec))) {
    cat("X bScore does not match Braak\n")
  } else if (verbose) {
    cat("OK Braak vs bScore\n")
  }
}


# Write a metadata data frame to a file
#
# This is a wrapper around `write.csv` that has some extra handling for values
# that contain commas and end of line characters (\n). Values with commas
# need to be escaped with quotes for a CSV file, and some data sets have them
# escaped already and some don't. \n characters are removed entirely, as quote
# escaping doesn't affect them.
#
# Arguments:
#   metadata - a data frame of harmonized metadata where rows are individuals
#              and columns are variables
#   filename - the base name of the original metadata file from Synapse (without
#              path information). This function automatically inserts
#              "_harmonized" just before ".csv" in the file name, and writes it
#              to data/output/
#
# Returns:
#   the new file name that was written, which should have "_harmonized" added
#   before ".csv"
write_metadata <- function(metadata, filename) {
  # Put quotes around values with commas
  for (column in colnames(metadata)) {
    if (is.character(metadata[, column])) {
      # Columns that contain commas and aren't already escaped with quotes
      commas <- grepl(",", metadata[, column]) & !grepl("\"", metadata[, column])
      metadata[commas, column] <- paste0("\"", metadata[commas, column], "\"")

      # Remove any "\n" characters
      metadata[, column] <- str_replace_all(metadata[, column], "\n", "")
    }
  }

  new_filename <- file.path(
    "data", "output", str_replace(filename, "\\.(csv|txt)", "_harmonized.csv")
  )
  write.csv(metadata, new_filename,
    row.names = FALSE, quote = FALSE
  )

  return(new_filename)
}


# Upload a file to Synapse
#
# This is a wrapper around `synStore` to shorten code slightly. It uploads a
# file to Synapse but only if the contents are different than what is currently
# on Synapse. This check is done because if the file is the same, synStore will
# erase any version comments in Synapse even if the file itself doesn't get a
# new version number. This function also makes sure that synStore doesn't erase
# any annotations on the file in Synapse.
#
# Arguments:
#   filename - the full path and name of the file to upload
#   folder_id - the Synapse ID of the folder on Synapse where the file should be
#     uploaded.
#   prov_used - (optional) a string or vector of Synapse IDs or URLs of files
#     that were used as input, for provenance.
#   prov_executed - (optional) a string or vector of Synapse IDs or URLs (like
#     Github links) pointing to scripts that were run to process the data, for
#     provenance.
#
# Returns:
#   a Synapse `File` object containing information about the uploaded file
synapse_upload <- function(filename, folder_id, prov_used = NULL, prov_executed = NULL) {
  syn_info <- synapse_get_info(filename, folder_id)

  if (!is.null(syn_info)) {
    md5 <- tools::md5sum(filename)

    # Don't actually update the file.
    if (md5 == syn_info$get("_file_handle")$contentMd5) {
      message(str_glue("\"{basename(filename)}\" matches the file on Synapse ",
                       "and will not be re-uploaded."))
      return(syn_info)
    }
  }

  syn_file <- File(filename, parent = folder_id)
  syn_file <- synStore(syn_file, forceVersion = FALSE, set_annotations = FALSE,
                       used = prov_used, executed = prov_executed)
  return(syn_file)
}


# Search a folder on Synapse to see if a file is already there. If it is,
# return info on the file without downloading it. Otherwise, return NULL.
synapse_get_info <- function(filename, folder_id) {
  id <- synFindEntityId(basename(filename), folder_id)
  if (is.null(id)) {
    return(NULL)
  }
  synGet(id, downloadFile = FALSE)
}


# Download a file from Synapse
#
# This is a wrapper around `synGet` to shorten code slightly. All downloads go
# into "data/downloads", and if a file with that name already exists, the old
# file is overwritten with the new one to avoid making multiple copies.
#
# Arguments:
#   syn_id - the Synapse ID of the file on Synapse to download
#   dataset_name - optional, the name of the dataset. Used only for printing a
#     warning if the file has a newer version on Synapse.
#
# Returns:
#   a Synapse `File` object containing information about the downloaded file
synapse_download <- function(syn_id, dataset_name = NULL) {
  check_new_version(syn_id, dataset_name)

  synGet(syn_id,
         downloadLocation = file.path("data", "downloads"),
         ifcollision = "overwrite.local")
}


# Check Synapse for new file versions
#
# All Synapse IDs used for this code have the file's version number included for
# reproducibility. This function checks to see if there are newer versions
# available on Synapse than what is specified in the code, and prints a message
# if that's the case.
#
# Arguments:
#   syn_id - a Synapse ID of the format "syn123" or "syn123.5", where the
#     optional number after the decimal is the file version on Synapse. If no
#     version is specified, a warning is printed stating that the latest version
#     of the file on Synapse will be used.
#   dataset_name - optional, the name of the dataset. Used only for printing a
#     warning if the file has a newer version on Synapse.
#
# Returns:
#   nothing
check_new_version <- function(syn_id, dataset_name = NULL) {
  # Separate into ID and version
  vals <- str_split_1(syn_id, pattern = "\\.")

  # If no version given, print a warning
  if (length(vals) != 2 || is.na(suppressWarnings(as.numeric(vals[2])))) {
    warning(
      str_glue(
        "WARNING: no valid version specified for '{syn_id}' ({dataset_name}). ",
        "The latest version will be used for harmonization."
      ),
      "\n"
    )
  } else {
    # Version number exists
    id <- vals[1]
    version <- vals[2]
    syn_file <- synGet(id, downloadFile = FALSE)

    # If a newer version exists on Synapse, print a warning
    if (syn_file$versionNumber != version) {
      warning(
        str_glue(
          "WARNING: there is a new version of {id} ({dataset_name}): ",
          "{version} => {syn_file$versionNumber}. Version {version} will be ",
          "used for harmonization."
        ),
        "\n"
      )
    }

    # No message printed if the file matches the version on Synapse
  }
}
