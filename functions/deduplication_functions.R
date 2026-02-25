# This file contains the functions necessary to fill in missing information
# in studies that share the same individual, but may have differing levels of
# missing data for that individual.

library(dplyr)
library(stringr)
library(purrr)

# De-duplicate metadata from different studies
#
# This function takes a list of metadata data frames, concatenates them, and
# then attempts to resolve cases where different studies have different data for
# the same individual.
#
# For each individual that exists in multiple studies, resolve duplicates:
#   1. For columns where some rows have NA and some have a unique non-NA value,
#      replace the NA value with that unique non-NA value.
#   2. For columns where some rows have "missing or unknown" and some have a
#      unique value other than that, replace "missing or unknown" with the
#      unique value.
#   3. For the ageDeath/PMI columns where rows have different numbers, report
#      the difference but leave the values as-is. If rows disagree only because
#      of precision, nothing is reported.
#   4. For cases where cohort values disagree, resolve as follows:
#       a) Replace "ROSMAP" with the cohort value from Diverse Cohorts /
#          AMP-AD 1.0,
#       b) Otherwise, report the difference but don't change any values
#   5. For disagreements in the "apoeGenotype" or "apoeStatus" columns, use the
#      "NPS-AD" value if it exists, otherwise print it as unresolvable
#
# When duplicated data is un-resolvable, either because it is a special case
# that is intentionally flagged or because this function doesn't have anything
# implemented to handle it, the following information is printed out:
#   "<individualID> <column name> [<list of values for this individual/column>]"
#
# Note: To shorten this function and make it more readable, some processing has
# been broken out into separate functions.
#
# Arguments:
#   df_list - a list of data frames, each of which must include an `individualID`
#       column as well as every column listed in `include_cols`
#   overlap_individuals - a data.frame containing the individualIDs from each
#       study, and a grouping ID (groupID) to indicate overlap between studies,
#       as returned by `find_overlaps()`.
#   spec - a `config` object describing the standardized values for each field,
#       as defined by this project's `config.yml` file
#   include_cols - a vector of column names to de-duplicate
#   exclude_cols - a vector of column names that should be excluded from
#       de-duplication even if they are listed in `include_cols`
#   verbose - if FALSE, only issues or un-resolvable data will be printed. If
#       TRUE, information on every column with duplicate values for each
#       individual will be reported even if the duplication is resolved or ignored.
#
# Returns:
#   a single data frame containing all rows from all data frames in `df_list`,
#   with column values de-duplicated where possible. The data frame will contain
#   all columns present in any data frame in `df_list`, and values will be `NA`
#   for rows that come from data frames without that column.
deduplicate_studies <- function(df_list,
                                overlap_individuals,
                                spec,
                                include_cols = spec$required_columns,
                                exclude_cols = c("individualID", "dataContributionGroup"),
                                verbose = TRUE) {
  # Make sure certain fields in each data frame are of the same type, then merge
  # in the overlap data frame to be able to compare IDs across studies
  df_list <- lapply(df_list, function(df_item) {
    df_item |>
      mutate(across(any_of(c("individualID", "apoeGenotype", "amyAny")),
                    as.character)) |>
      # Adds "groupID" which ties individuals together across studies
      merge(overlap_individuals)
  })

  meta_all <- list_rbind(df_list)
  include_cols <- setdiff(include_cols, exclude_cols)

  # Find the IDs that exist in multiple studies
  dupe_ids <- overlap_individuals |>
    group_by(groupID) |>
    count() |>
    subset(n > 1) |>
    pull(groupID)

  # Resolve duplicated data for each individual as best as possible
  for (group_id in dupe_ids) {
    # This will be altered to resolve duplication, and will get added back to the
    # meta_all data frame
    meta_tmp <- subset(meta_all, groupID == group_id)

    for (col_name in include_cols) {
      unique_vals <- unique(meta_tmp[, col_name])

      if (length(unique_vals) > 1) {
        # For reporting un-resolved mismatches
        report_string <- paste(
          group_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n"
        )

        # This prints out a mismatch even if it's resolved below. If verbose is
        # FALSE, only unresolvable mismatches will be printed.
        if (verbose) {
          cat(report_string)
        }

        # Remove NA, "", and "missing or unknown" values and see what's left
        leftover <- setdiff(unique_vals, c(spec$missing, "")) |>
          na.omit()

        # If one unique value left, replace all values with that value
        if (length(leftover) == 1) {
          meta_tmp[, col_name] <- leftover
        } else if (length(leftover) == 0) {
          # Nothing left, use "missing or unknown" if it's there, otherwise set
          # to NA.
          if (any(unique_vals == spec$missing)) {
            meta_tmp[, col_name] <- spec$missing
          } else {
            meta_tmp[, col_name] <- NA
          }
        } else {
          # More than 1 value left after omission of NA and "missing or unknown".
          # Do some column-specific handling of duplicates
          if (col_name %in% c("apoeGenotype", "apoe4Status") &&
              "NPS-AD" %in% meta_tmp$study &&
              meta_tmp[meta_tmp$study == "NPS-AD", col_name] != spec$missing) {
            # Use the NPS-AD value for apoe genotype / status where it disagrees
            # with other data sets
            meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
          } else if (col_name %in% c("ageDeath", "PMI")) {
            meta_tmp <- deduplicate_ageDeath_pmi(meta_tmp, leftover, col_name, report_string)
          } else if (col_name == "cohort") {
            meta_tmp <- deduplicate_cohort(meta_tmp, leftover, col_name, spec, report_string)
          } else {
            # Column is something else that we don't have specific handling for,
            # so we report it but don't try to resolve duplication.
            cat(report_string)
          }
        }
      }
    }

    # Replace original rows with de-duplicated data
    rows_replace <- meta_all$groupID == group_id
    meta_all[rows_replace, ] <- meta_tmp
  }

  return(meta_all)
}


# Age- and PMI-specific handling for de-duplication
#
# This function is used inside `deduplicate_studies` to resolve duplication of
# age and PMI data for a single individual. If this function is called, then
# there are at least 2 distinct, non-NA values in the `ageDeath` or `PMI`
# column that are assigned to this individual. This function checks whether the
# numbers differ only by precision (e.g. 3.5 vs 3.547), or if the numbers are
# completely different (e.g. 4 vs 10). The former case is ignored, and the
# latter case will be reported to the console. Currently, no modification is
# done to the values themselves, as NPS-AD wishes to keep all of their values
# as-is even where they differ from Diverse Cohorts / AMP-AD 1.0.
#
# Arguments:
#   meta_tmp - a data frame with 2 or more rows, where all rows have the same
#     individual ID and there are 2 or more distinct values in the `ageDeath`
#     and/or `PMI` column
#   leftover - a vector of unique values from <col_name> for this individual,
#     which has had `NA` values removed
#   col_name - the name of the column being handled, either "ageDeath" or "PMI"
#   report_string - the string that gets printed out if there is a real
#     difference between values that is not due to precision. The string
#     contains the `individualID`, column name, and unique values in the column.
#
# Returns:
#   meta_tmp unaltered (giving us the option to alter it in a future update)
deduplicate_ageDeath_pmi <- function(meta_tmp, leftover, col_name, report_string) {
  # By this point there is more than one unique, non-NA age value in `leftover`.
  # Check if all values are roughly equal to make sure there isn't an actual
  # mis-match. There may be more than 2 values so this is the quickest way to
  # check.
  num_vals <- suppressWarnings(as.numeric(leftover)) |>
    na.omit()

  equivalent <- sapply(num_vals, all.equal, num_vals[1], tolerance = 1e-3)

  # If there's a real mismatch, report it. If a value is missing and there are
  # multiple studies with equivalent but not identical values, first try using
  # the NPS-AD value if it exists, and then the Diverse Cohorts value. This may
  # need to be updated if more overlapping studies are added.
  # Note that `all.equal` returns a string with the difference between two
  # numbers if they are not equal, rather than FALSE, so we have to check for !=
  # TRUE instead of == FALSE or !equivalent.
  if (any(equivalent != TRUE)) {
    cat(report_string)
  } else if (length(leftover) != ncol(meta_tmp)) {
    # No real mismatch but at least one value was NA and there are multiple
    # close-enough values. Use NPS-AD value first if it exists, then Diverse
    # Cohorts if it exists. Otherwise print.
    na_vals <- which(is.na(meta_tmp[, col_name]))

    if ("NPS-AD" %in% meta_tmp$study) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
    } else if ("AMP-AD_DiverseCohorts" %in% meta_tmp$study) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "AMP-AD_DiverseCohorts", col_name]
    } else {
      cat("Unresolved NA fill: ", report_string)
    }
  }

  # Otherwise we leave the ages as-is even if there's difference in precision

  return(meta_tmp)
}


# Cohort-specific handling for de-duplication
#
# This function is used inside `deduplicate_studies` to resolve duplication of
# cohort values for a single individual. If this function is called, then
# there are at least 2 distinct, non-NA values in the `cohort` column that are
# assigned to this individual. This function handles one special case and
# defaults to just reporting the duplication if the special case doesn't apply:
#   * NPS-AD reports the `cohort` of all ROSMAP samples as "ROSMAP" rather than
#     identifying them as "ROS" or "MAP" separately. All of these samples exist
#     in Diverse Cohorts or ROSMAP 1.0 metadata, which has the correct
#     separation into "ROS" and "MAP", so we replace NPS-AD's ROSMAP `cohort`
#     values with the correct values from DC/ROSMAP 1.0.
#
# Arguments:
#   meta_tmp - a data frame with 2 or more rows, where all rows have the same
#     individual ID and there are 2 or more distinct values in the `cohort`
#     column
#   leftover - a vector of unique values from <col_name> for this individual,
#     which has had `NA` values removed
#   col_name - the name of the column being handled ("cohort")
#   report_string - the string that gets printed out if there is a real
#     difference between values that is not due to precision. The string
#     contains the `individualID`, column name, and unique values in the column.
#
# Returns:
#   meta_tmp, which will have some `cohort` values replaced if the original
#   value was "ROSMAP". Other values and columns are left as-is.
deduplicate_cohort <- function(meta_tmp, leftover, col_name, spec, report_string) {
  # If there is more than one left over value and the values don't meet the
  # special case below, report it but don't try to resolve duplication.

  # Special case: NPS-AD reports cohort on some samples as "ROSMAP", which needs
  # to instead use the cohort value from Diverse Cohorts or ROSMAP 1.0 data.
  if ("ROSMAP" %in% leftover) {
    cohort_val <- setdiff(leftover, "ROSMAP")

    # If there is still more than one value for cohort, or no remaining values,
    # report it but don't resolve the de-duplication
    if (length(cohort_val) != 1) {
      cat(report_string)
    } else {
      # Otherwise resolve duplication
      meta_tmp[, col_name] <- cohort_val
    }
  } else {
    # Special case doesn't apply, report the mismatch but don't resolve it.
    cat(report_string)
  }

  return(meta_tmp)
}


# Fill missing AMP-AD 1.0 values in Diverse Cohorts
#
# The Diverse Cohorts metadata has an `individualID_AMPAD_1.0` column that is
# supposed to have the corresponding AMP-AD 1.0 ID if that sample exists in 1.0
# data. However this field has a lot of NAs for samples that do exist in 1.0
# data, whether intentional or not. This function fills in the matching ID from
# 1.0 data in these cases.
#
# Arguments:
#   meta_all - a data.frame of all harmonized study data, as returned by `deduplicate_studies()`.
#   overlap_individuals - a data.frame containing the individualIDs from each
#     study, and a grouping ID (groupID) to indicate overlap between studies,
#     as returned by `find_overlaps()`.
#   studies - the "studies" object from config.yml
#
# Returns:
#   meta_all but with appropriate `individualID_AMPAD_1.0` values filled in
fill_missing_ampad1.0_ids <- function(meta_all, overlap_individuals, studies) {
  dc <- subset(meta_all, study == studies$diverse_cohorts$name) |>
    select(individualID, cohort, individualID_AMPAD_1.0, study) |>
    merge(overlap_individuals)

  matches_1.0 <- subset(overlap_individuals, groupID %in% dc$groupID &
                          study %in% c(studies$mayo$name, studies$msbb$name,
                                       studies$rosmap$name)) |>
    dplyr::rename(individualID_1.0_tmp = individualID) |>
    select(-study)

  # There shouldn't be any duplicated unique IDs since the three 1.0 studies
  # don't have any overlap
  stopifnot(length(unique(matches_1.0$groupID)) == nrow(matches_1.0))

  missing_vals <- is.na(dc$individualID_AMPAD_1.0) &
    dc$groupID %in% matches_1.0$groupID

  cat(str_glue("Filling {sum(missing_vals)} missing AMPAD-1.0 IDs in Diverse Cohorts"), "\n")

  # Merge in the 1.0 individual IDs from matches_1.0 as a column so we can do
  # direct replacement
  dc <- merge(dc, matches_1.0, all = TRUE) |>
    mutate(individualID_AMPAD_1.0 = ifelse(is.na(individualID_AMPAD_1.0),
                                           individualID_1.0_tmp,
                                           individualID_AMPAD_1.0)) |>
    select(-individualID_1.0_tmp)

  col_order <- colnames(meta_all)

  meta_all |>
    # Replace column
    select(-individualID_AMPAD_1.0) |>
    merge(dc, all = TRUE, sort = FALSE) |>
    # Restore original column order
    select(all_of(col_order))
}
