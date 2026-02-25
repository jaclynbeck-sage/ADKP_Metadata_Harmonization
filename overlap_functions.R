library("igraph")

# Find overlapping individuals between studies
#
# This is the "main" overlap finding function. The overlap-finding process was
# broken out into several functions for readability, and this function is the
# one that is called to kick off the process.
#
# The main process is:
#   1. Create a copy of the harmonized metadata and alter the format of some
#      studies' individualIDs to make them comparable across studies.
#   2. Get all (modified) IDs across all studies and find the ones that exist in
#      more than one study.
#   3. For every pair of studies that share an ID, compare values for that ID in
#      each overlapping column between the pair. Values are marked as "match" if
#      they are equivalent, "missing" if one or both values are missing, or
#      "mismatch" if the values do not match AND neither of them is missing.
#   4. For every pair of studies that share an ID, count up the number of
#      matches and mismatches for that ID across all shared columns between
#      the pair.
#   5. Determine which studies with the same ID actually represent the same individual:
#       * A pair of studies definitely share the same individual if:
#         * There are no mismatches between *non-missing* values in overlapping columns
#         * Less than 50% of the data is missing one or both values in a given column
#       * Some special case individuals from NPS-AD and ROSMAP that do not meet
#         these criteria were verified as matching by hand and added to the list
#         of true overlaps.
#   6. Create a data frame that maps overlapping studies' individualIDs (which
#      may have different formatting) to the same "individual". This is done by
#      creating a "groupID" that represents a single individual, and mapping
#      that groupID to the real individualID in each study with that individual.
#      IDs from a single study that have no overlap will have a groupID unique
#      to that one study.
#
# The output data frame will look like this:
#
#   groupID   study   individualID
#   123;1     A       123
#   123;1     B       AMPAD_MSSM_000123
#   123;2     C       123
#   456;1     C       456
#
# In this example, studies A and B share an individual (same groupID), but have
# differently-formatted individualIDs. Study C has an individual with the same
# ID (123) as studies A and B, but it's *not* the same individual, so study C
# has a different groupID for that individual than studies A and B. Study C also
# has an individual (456) that doesn't exist in the other studies so it has a
# unique groupID too.
#
# Arguments:
#   df_list - a named list of data frames, which contain the harmonized metadata
#     for each study. The name of each list item should be the study's name as
#     defined in config.yml
#   studies - the "studies" object from config.yml
#   spec - the "columns" object from config.yml
#
# Returns:
#   a data frame with columns "groupID", "study", and "individualID", where
#   individualID is the ID of the individual as it appears in the study, and
#   groupID is a value that identifies groups of studies/individualIDs that
#   represent the same individual. An example is defined above in the function
#   description.
#
find_overlaps <- function(df_list, studies, spec) {
  study_names <- names(df_list)

  # Step 1: Modify all of the study data to use consistent ID formats and handle
  # a few special cases.
  #   * MSBB IDs need to have "AMPAD_MSSM_00.." removed to match Diverse Cohorts IDs.
  #   * Samples in Diverse Cohorts that have dataContributionGroup = Emory, but
  #     these samples should pair up with MSBB samples. We change the
  #     dataContributionGroup to MSSM so these samples match between Diverse
  #     Cohorts and MSBB 1.0.
  #   * Special case: MC_snRNA uses "Control" in its diagnosis field instead of
  #     "control", which is used in every other Mayo study. We change the values
  #     to lower-case to make comparison with other Mayo studies easier.
  # These changes apply to the overlap finding functions only and do not end up
  # in the final harmonized files.
  mod_data <- lapply(df_list, function(study_df) {
    study_name <- unique(study_df$study)

    data <- study_df |>
      mutate(orig_individualID = as.character(individualID),
             individualID = msbb_ids_to_divco(individualID, dataContributionGroup),
             # Round all numeric columns to 2 digits for comparison
             across(where(is.numeric), ~round(.x, digits = 2)),
             # Round numeric ageDeath values too
             ageDeath = ifelse(
               ageDeath == spec$ageDeath$over90,
               ageDeath,
               suppressWarnings(as.numeric(ageDeath)) |> round(digits = 2) |> as.character()
             ),
             # Change Emory dataContributionGroup values to MSSM
             dataContributionGroup = ifelse(dataContributionGroup == spec$dataContributionGroup$emory,
                                            spec$dataContributionGroup$mssm,
                                            dataContributionGroup),
             # Used for group finding later on
             node_id = paste0(study_name, ";", individualID)) |>
      select(-filename)

    # Special case: MC_snRNA "Control" -> "control"
    if (study_name == studies$mc_snrna$name) {
      data$diagnosis[data$diagnosis == "Control"] <- "control"
    }

    return(data)
  })

  # Step 2: Get all IDs and which study they come from
  all_ids <- lapply(mod_data, select, individualID, study, orig_individualID, node_id) |>
    list_rbind()

  # Step 3: Create a data frame containing info on whether values for the same
  # ID in multiple studies match or not
  match_df <- create_match_df(all_ids, mod_data, studies, spec)

  # Steps 4-6: Determine which overlapping IDs represent the same individual
  match_individuals(match_df, all_ids, studies)
}


# Find columns present in both data sets
#
# This function finds columns that appear in both data sets for the purpose of
# calculating equivalence of two rows of data. This requires excluding certain
# columns to avoid skewing the comparison numbers. For example, IDs and study
# names will always mismatch, and derived value matches/mismatches are redundant
# with match/mismatches in the columns they derive from.
#
# This function also handles several special cases where a column is known to be
# completely empty and shouldn't count against the missing data stats for that
# data set to avoid skewing the numbers:
#   * ROSMAP 1.0 data will never have any amyThal values
#   * Mayo 1.0 data will never have any amyCerad values
#   * NPS-AD race and isHispanic values were all set to missing during harmonization
#
# Arguments:
#   pair - a character vector with the names of the two studies being compared,
#     which is necessary to make sure study-specific special cases are handled
#   data1 - a data frame of metadata from the first study being compared
#   data2 - a data frame of metadata from the second study being compared
#
# Returns:
#   a character vector of column names that overlap between data sets, excluding
#   any special cases mentioned above
#
find_matching_columns <- function(pair, data1, data2) {
  same_cols <- intersect(colnames(data1), colnames(data2)) |>
    # Exclude IDs and derived columns
    setdiff(c("individualID", "study", "orig_individualID", "node_id",
              "species", "apoe4Status", "ADoutcome", "amyA", "amyAny",
              "bScore"))

  # Special case: ROSMAP 1.0 data will not have any amyThal values so this
  # shouldn't count against matching
  if (any(pair == studies$rosmap$name)) {
    same_cols <- setdiff(same_cols, "amyThal")
  }

  # Special case: Mayo 1.0 data will not have any amyCerad values so this
  # shouldn't count against matching
  if (any(pair == studies$mayo$name)) {
    same_cols <- setdiff(same_cols, "amyCerad")
  }

  # Special case: all NPS-AD race and isHispanic values are all missing and
  # shouldn't count against matching.
  if (any(pair == studies$nps_ad$name)) {
    same_cols <- setdiff(same_cols, c("race", "isHispanic"))
  }

  return(same_cols)
}


# Create a data frame to compare data for the same individualID
#
# This function looks at each pair of studies that contain the same individualID
# and creates a data frame denoting whether, for a given column, the two values
# for that individual match, have one missing value, have both values missing,
# or have a mismatch. The set of columns being compared is the set of
# overlapping columns between the two studies (all harmonized columns plus
# potentially some un-harmonized columns), minus any special cases described in
# the find_matching_columns function.
#
# Arguments:
#   all_ids - a data frame containing all IDs across all the studies
#   study_dfs - a named list of all data frames containing harmonized metadata
#     for all studies. The name of each list item should be the study's name as
#     defined in config.yml
#   studies - the "studies" object from config.yml
#   spec - the "columns" object from config.yml
#
# Returns:
#   a data frame with one row per id + study pair + column name, containing
#   which studies were compared, the column name, each study's value for that
#   column/individual, and whether the values match, are missing, or mismatch.
#   There is more information in this data frame than strictly necessary to
#   determine overlap, in order to make debugging and hand comparison easier.
#
create_match_df <- function(all_ids, study_dfs, studies, spec) {
  # Find IDs appearing in more than one dataset
  matched_ids <- table(all_ids$individualID)
  matched_ids <- names(matched_ids)[matched_ids > 1]

  # Studies that have overlapping IDs but we already know don't have overlapping
  # samples. Pairs are sorted alphabetically so setdiff() works in the lapply
  # below.
  bad_pairs <- list(
    c(studies$mayo$name, studies$nps_ad$name) |> sort(), # Mayo / NPS-AD
    c(studies$mayo$name, studies$msbb$name) |> sort(), # Mayo / MSBB
    c(studies$mc_snrna$name, studies$msbb$name) |> sort() # MC_snRNA / MSBB
  )

  # For each ID that appears in at least two studies, create a data frame with
  # the combined info for all pairs of studies it appears in
  match_df <- lapply(matched_ids, function(id) {
    studies_check <- subset(all_ids, individualID == id) |> pull(study)

    # All possible pairs from studies_check, removing bad pairs we already
    # know about
    study_pairs <- combn(studies_check, m = 2, simplify = FALSE) |>
      setdiff(bad_pairs)

    # Create one data frame per pair and rbind all data frames together
    lapply(study_pairs, function(pair) {
      data1 <- subset(study_dfs[[pair[1]]], individualID == id)
      data2 <- subset(study_dfs[[pair[2]]], individualID == id)

      same_cols <- find_matching_columns(pair, data1, data2)

      data.frame(id = id,
                 study1 = pair[1],
                 study2 = pair[2],
                 node_id1 = data1$node_id,
                 node_id2 = data2$node_id,
                 column = same_cols,
                 val1 = as.character(t(data1[, same_cols])),
                 val2 = as.character(t(data2[, same_cols])))
    }) |>
      list_rbind()
  })

  # Bind all data frames for all IDs together and check for value matches in
  # each row
  match_df <- list_rbind(match_df) |>
    mutate(
      # Convert all "missing or unknown" values to NA for easier missingness checks
      val1 = na_if(val1, spec$missing),
      val2 = na_if(val2, spec$missing),
      # Equivalence or missing values
      match = val1 == val2,
      both_missing = is.na(val1) & is.na(val2),
      one_missing = (is.na(val1) | is.na(val2)) & !both_missing,
      # Fix NA match values (which happen when either val1 or val2 are NA)
      match = ifelse(is.na(match), FALSE, match)
    )

  # Special case: NPS-AD cohort may be "ROSMAP" instead of "ROS" or "MAP"
  # but we consider that a match. The "if" statement is TRUE if a) NPS-AD is
  # one of the studies, and b) one of the cohort values is "ROSMAP"
  nps_cohort <- (match_df$study1 == studies$nps_ad$name |
                   match_df$study2 == studies$nps_ad$name) &
    match_df$column == "cohort" &
    (match_df$val1 == "ROSMAP" | match_df$val2 == "ROSMAP")

  match_df$match[nps_cohort] <- TRUE

  # Anything that isn't a match or missing data is a mismatch
  match_df$mismatch <- !match_df$match & !match_df$one_missing & !match_df$both_missing
  return(match_df)
}


# Steps 4-6: Determine if individualIDs that match across studies really
# represent the same individual
#
# Arguments:
#   match_df - a data frame as output by `create_match_df`, with one row per
#     combination of ID, study pair, and column, describing whether the values
#     match or not
#   all_ids - a data frame containing columns for (modified) individualID, study
#     name, original individualID, and a "node ID" that is used for clique
#     finding
#   studies - the "studies" object from config.yml
#
# Returns:
#   a data frame with columns "groupID", "study", and "individualID", where
#   individualID is the ID of the individual as it appears in the study, and
#   groupID is a value that identifies groups of studies/individualIDs that
#   represent the same individual.
match_individuals <- function(match_df, all_ids, studies) {
  # Step 4: Count up all the matches/missing data/mismatches across all columns
  # for each individual in each pair of studies. Some extra information is
  # included in this data frame to make debugging and hand comparison easier.
  stats <- match_df |>
    group_by(id, study1, study2, node_id1, node_id2) |>
    summarize(
      n_match = sum(match),
      n_one_missing = sum(one_missing),
      n_both_missing = sum(both_missing),
      n_mismatch = sum(mismatch),
      # Don't count columns where both studies are missing data in the
      # percentage calculation for matches
      pct_match = n_match / (length(match) - n_both_missing),
      fields_match = paste(column[match & !both_missing], collapse = " / "),
      fields_one_missing = paste(column[one_missing], collapse = " / "),
      fields_both_missing = paste(column[both_missing], collapse = " / "),
      fields_mismatch = paste(column[mismatch], collapse = " / "),
      .groups = "drop"
    )

  # Step 5: Two studies definitely refer to the same individual if they have a
  # match percent higher than 50% and no mismatches. This indicates that up to
  # 50% of the data can be missing, but the non-missing data matches 100%.
  same_samples <- subset(stats, pct_match >= 0.5 & n_mismatch == 0)

  # Special case: There is 1 ROSMAP sample that is missing more than 50% of its
  # data but has been verified by hand as a match with NPS-AD. We run this
  # generic check instead of hard-coding the sample ID so that we are aware if
  # this happens with any other new studies.
  missing_data_samples <- subset(stats, pct_match < 0.5 & n_mismatch == 0 &
                                   grepl("cohort", fields_match))

  # Check: There should only be 1 row, paired between ROSMAP and NPS-AD. If more
  # rows show up after adding a new study, this will throw an error so we can
  # verify by hand if the new individual matches, and can alter this statement
  # so it passes again.
  nps_name <- studies$nps_ad$name
  stopifnot(
    nrow(missing_data_samples) == 1 &&
      all(c(nps_name, studies$rosmap$name) %in%
            c(missing_data_samples$study1, missing_data_samples$study2))
  )

  # Special case: Until de-duplication, NPS-AD data might have up to 1 mismatch
  # with Diverse Cohorts or MSBB data but will still have > 50% matches. The
  # samples in this subset have been verified as matching by hand.
  nps_samples <- subset(stats,
                        (study1 == nps_name | study2 == nps_name) &
                          pct_match >= 0.5 & n_mismatch == 1)

  # Check: There should only be 8 samples, and the only studies present should
  # be Diverse Cohorts, MSBB, and NPS-AD. If this changes after adding a new
  # study, the new samples need to be verified by hand and added to the check
  # so it passes.
  stopifnot(nrow(nps_samples) == 8 &&
              all(c(nps_samples$study1, nps_samples$study2) %in%
                    c(studies$diverse_cohorts$name, studies$msbb$name, nps_name)))

  # Add the hand-verified missing data samples and NPS samples to the list of
  # matches
  same_samples <- rbind(same_samples, missing_data_samples) |>
    rbind(nps_samples)

  # Anything with more than 1 mismatch is not a matching sample. This happens
  # because Mayo and numerical MSBB IDs overlap in Diverse Cohorts.
  no_match <- subset(stats, n_mismatch > 1)

  # Make sure we've covered all rows in stats with either same_samples or
  # no_match, and make sure there are no overlapping rows between same_samples
  # and no_match. If this check fails, then we have missed a special case and
  # need to investigate what that is.
  stopifnot(nrow(same_samples) + nrow(no_match) == nrow(stats) &&
              nrow(distinct(rbind(same_samples, no_match))) == nrow(stats))

  # Step 6: Map studies with the same individual to the same group.

  # Use the same_samples df to find all overlapping groups of studies for each
  # ID. This can be done efficiently by treating this as an undirected graph and
  # finding all max cliques. Each node in the graph is "<study>;<id>" with edges
  # between each pair that have the same individual. max_cliques will find the
  # largest group of studies for each ID that are all connected to each other.
  # For example, if an individualID is present in pairs of studies represented
  # by A-B, B-C, A-C, and D-E, studies A, B, and C are a clique and contain the
  # same individual, while studies D and E are a separate clique representing a
  # different individual that happens to have the same ID.
  graph <- same_samples |>
    select(node_id1, node_id2) |>
    igraph::graph_from_data_frame(directed = FALSE)

  groups <- igraph::max_cliques(graph)

  # Put the cliques list back into data frame form, with one group ID per
  # clique
  groups <- lapply(1:length(groups), function(group_id) {
    data.frame(group_id = as.character(group_id),
               node_id = igraph::as_ids(groups[[group_id]]))
  }) |>
    list_rbind()

  # Merge the grouping information back into the data frame with all IDs,
  # matching on node_id. IDs in all_ids that don't exist in the group data must
  # only exist in a single study, and end up with an NA group_id.
  group_df <- merge(all_ids, groups, all = TRUE)

  # Add a "none" group for IDs that don't belong to a group, but make each value
  # unique so non-matched entries with the same ID don't get grouped together
  na_vals <- which(is.na(group_df$group_id))
  group_df$group_id[na_vals] <- paste0("none.", 1:length(na_vals))

  # Assign a unique ID for each *group* of overlapping individualIDs, so that
  # non-overlapping studies with the same ID can be distinguished. Then, replace
  # the individualID values with the original, un-modified individualIDs from
  # each study. The group_id column contains numbers when a group exists, and we
  # convert this to a string based on the shared individualID to avoid issues
  # with indexing or numeric / string collisions. We use
  # as.numeric(factor(group_id)) for cosmetic reasons, so that no matter what
  # values exist in group_id for a given individualID, they are assigned numbers
  # that start at 1.
  group_df <- group_df |>
    group_by(individualID) |>
    mutate(groupID = paste0(individualID, ";", as.numeric(factor(group_id)))) |>
    ungroup() |>
    select(groupID, study, orig_individualID) |>
    dplyr::rename(individualID = orig_individualID)

  return(group_df)
}
