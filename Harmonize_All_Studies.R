library(purrr)
library(synapser)
library(dplyr)
library(stringr)
library(readxl)

source(file.path("functions", "deduplication_functions.R"))
source(file.path("functions", "harmonize_functions.R"))
source(file.path("functions", "overlap_functions.R"))
source(file.path("functions", "util_functions.R"))

studies <- config::get("studies")
spec <- config::get("columns")

# TRUE will print column names of metadata + unique values for each dataset.
# FALSE will only print the status of the harmonized metadata for each dataset.
verbose <- FALSE

synLogin()

# Pre-initialize the list where harmonized data will go
df_list <- vector(mode = "list", length = length(studies))
names(df_list) <- sapply(studies, "[[", "name")


# Diverse Cohorts --------------------------------------------------------------

diverse_cohorts <- studies$diverse_cohorts

meta_new <- harmonize(diverse_cohorts, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", diverse_cohorts$name, "\n")
validate_values(meta_new, spec)

df_list[[diverse_cohorts$name]] <- meta_new


# MayoRNAseq -------------------------------------------------------------------

mayo <- studies$mayo

# Extra metadata - the biospecimen file from Synapse
biospec_file <- synapse_download(mayo$extra_metadata, mayo$name)
biospecimen <- read.csv(biospec_file$path)

meta_new <- harmonize(mayo, spec, extra_metadata = biospecimen)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", mayo$name, "\n")
validate_values(meta_new, spec)

df_list[[mayo$name]] <- meta_new


# MC-BrAD ----------------------------------------------------------------------

mc_brad <- studies$mc_brad

meta_new <- harmonize(mc_brad, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", mc_brad$name, "\n")
validate_values(meta_new, spec)

df_list[[mc_brad$name]] <- meta_new


# MC_snRNA ---------------------------------------------------------------------

mc_snrna <- studies$mc_snrna

meta_new <- harmonize(mc_snrna, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", mc_snrna$name, "\n")
validate_values(meta_new, spec)

df_list[[mc_snrna$name]] <- meta_new


# MCMPS ------------------------------------------------------------------------

mcmps <- studies$mcmps

meta_new <- harmonize(mcmps, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", mcmps$name, "\n")
validate_values(meta_new, spec)

df_list[[mcmps$name]] <- meta_new


# MSBB -------------------------------------------------------------------------

msbb <- studies$msbb

meta_new <- harmonize(msbb, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", msbb$name, "\n")
validate_values(meta_new, spec)

df_list[[msbb$name]] <- meta_new


# NPS-AD -----------------------------------------------------------------------

nps <- studies$nps_ad

# Extra metadata -- neuropathology file from Synapse
neuro_file <- synapse_download(nps$extra_metadata, nps$name)
neuropath <- read.csv(neuro_file$path) |>
  dplyr::rename(individualID = IndividualID)

meta_new <- harmonize(nps, spec, extra_metadata = neuropath)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", nps$name, "\n")
validate_values(meta_new, spec) # Cohort will fail until after harmonization

df_list[[nps$name]] <- meta_new


# ROSMAP -----------------------------------------------------------------------

rosmap <- studies$rosmap

meta_new <- harmonize(rosmap, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", rosmap$name, "\n")
validate_values(meta_new, spec)

df_list[[rosmap$name]] <- meta_new


# SEA-AD -----------------------------------------------------------------------

# The version of SEA-AD that is on Synapse is missing Hispanic/Latino
# information that is present in the version released by the Allen Institute on
# brain-map.org. We use the version on Synapse but pull in the missing
# information from the Allen Institute version.

sea_ad <- studies$sea_ad

# Extra metadata -- Excel file from brain-map.org
sea_ad_file <- file.path("data", "downloads", basename(sea_ad$extra_metadata))
download.file(sea_ad$extra_metadata, destfile = sea_ad_file)

meta_sea_ad <- read_xlsx(sea_ad_file)

meta_new <- harmonize(sea_ad, spec, extra_metadata = meta_sea_ad)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", sea_ad$name, "\n")
validate_values(meta_new, spec)

df_list[[sea_ad$name]] <- meta_new


# SMIB-AD ----------------------------------------------------------------------

smib_ad <- studies$smib_ad

meta_new <- harmonize(smib_ad, spec)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", smib_ad$name, "\n")
validate_values(meta_new, spec)

df_list[[smib_ad$name]] <- meta_new


# Merge all files into one data frame ------------------------------------------

overlap_individuals <- find_overlaps(df_list, studies, spec)
overlap_file <- file.path("data", "output", "study_overlap_individualIDs.csv")
write.csv(arrange(overlap_individuals, groupID, study),
          overlap_file, row.names = FALSE, quote = FALSE)

meta_all <- deduplicate_studies(df_list, overlap_individuals, spec, verbose = FALSE)

meta_all <- meta_all |>
  fill_missing_ampad1.0_ids(overlap_individuals, studies) |> # For Diverse Cohorts
  updateADoutcome(spec) # For Diverse Cohorts

if (verbose) {
  cat("\nAll studies (de-duplicated)\n")
  print_qc(meta_all, spec)
}

cat("\nAll studies (de-duplicated)\n")
validate_values(meta_all, spec)


# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
new_files <- lapply(df_list, function(meta_old) {
  study_name <- unique(meta_old$study)

  meta_new <- subset(meta_all, study == study_name) |>
    select(all_of(colnames(meta_old))) |>
    # Remove the filename and study columns we added
    select(-filename, -study)

  # The order of individuals in a data set can get changed across the
  # deduplication functions so we make sure the individuals are back in the
  # original order. This makes it easier to compare across revisions.
  stopifnot(length(intersect(meta_old$individualID, meta_new$individualID)) ==
              length(unique(meta_new$individualID)))
  rownames(meta_new) <- meta_new$individualID
  meta_new <- meta_new[meta_old$individualID, ]

  return(list(study = study_name,
              filename = write_metadata(meta_new, unique(meta_old$filename))))
})

# Upload to ADKP metadata space

upload_id <- config::get("upload_synID")
github <- config::get("github_url")

# Metadata files
for (item in new_files) {
  # Get provenance
  which_study <- which(sapply(studies, "[[", "name") == item$study)
  study_obj <- studies[[which_study]]

  prov_used <- study_obj$metadata_id

  if ("extra_metadata" %in% names(study_obj)) {
    prov_used <- c(prov_used, study_obj$extra_metadata)
  }

  synapse_upload(item$filename, upload_id, prov_used, prov_executed = github)
}

# Study overlap file
prov_used <- sapply(studies, "[[", "metadata_id")
synapse_upload(overlap_file, upload_id, prov_used, prov_executed = github)
