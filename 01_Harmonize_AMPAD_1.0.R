# This script harmonizes and de-duplicates AMP-AD 1.0 metadata (MayoRNAseq,
# MSBB, and ROSMAP) with Diverse Cohorts and NPS-AD metadata, filling in missing
# information where there is sample overlap between the five files. The
# de-duplicated data is then used as the standard to check quality of other
# overlapping ADKP studies and fill in missing information in those studies
# as well.
#
# TODO add provenance

library(purrr)
library(synapser)
library(dplyr)
library(stringr)

source("util_functions.R")
source("harmonize_function.R")

studies <- config::get("studies")
spec <- config::get("columns")

# TRUE will print column names of metadata + unique values for each dataset.
# FALSE will only print the status of the harmonized metadata for each dataset.
verbose <- FALSE

synLogin()

df_list <- list()


# MayoRNAseq -------------------------------------------------------------------

mayo <- studies$mayo
meta_file <- synapse_download(mayo$metadata_id, mayo$name)
meta <- read.csv(meta_file$path)

biospec_file <- synapse_download(mayo$extra_metadata, mayo$name)
biospecimen <- read.csv(biospec_file$path)

if (verbose) {
  cat("\n", mayo$name, "\n")
  print(colnames(meta))
  print_qc(meta, spec,
           column_renames = list(
             PMI = "pmi",
             isHispanic = "ethnicity",
             amyCerad = "CERAD",
             amyThal = "Thal"
           )
  )
}

meta_new <- harmonize(mayo$name, meta, spec, extra_metadata = biospecimen) |>
  mutate(filename = meta_file$name)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", mayo$name, "\n")
validate_values(meta_new, spec)

df_list[["MayoRNAseq"]] <- meta_new


# MSBB -------------------------------------------------------------------------

msbb <- studies$msbb
meta_file <- synapse_download(msbb$metadata_id, msbb$name)
meta <- read.csv(meta_file$path)

if (verbose) {
  cat("\n", msbb$name, "\n")
  print(colnames(meta))
  print_qc(meta, spec,
           column_renames = list(
             PMI = "pmi",
             isHispanic = "ethnicity",
             amyCerad = "CERAD"
           )
  )
}

meta_new <- harmonize(msbb$name, meta, spec) |>
  mutate(filename = meta_file$name)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", msbb$name, "\n")
validate_values(meta_new, spec)

df_list[["MSBB"]] <- meta_new


# ROSMAP -----------------------------------------------------------------------

rosmap <- studies$rosmap
meta_file <- synapse_download(rosmap$metadata_id, rosmap$name)
meta <- read.csv(meta_file$path)

if (verbose) {
  cat("\n", rosmap$name, "\n")
  print(colnames(meta))
  print_qc(meta, spec,
           column_renames = list(
             ageDeath = "age_death",
             PMI = "pmi",
             sex = "msex",
             isHispanic = "spanish",
             apoeGenotype = "apoe_genotype",
             Braak = "braaksc",
             amyCerad = "ceradsc"
           )
  )
}

meta_new <- harmonize(rosmap$name, meta, spec) |>
  mutate(filename = meta_file$name)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", rosmap$name, "\n")
validate_values(meta_new, spec)

df_list[["ROSMAP"]] <- meta_new


# Diverse Cohorts --------------------------------------------------------------

diverse_cohorts <- studies$diverse_cohorts
meta_file <- synapse_download(diverse_cohorts$metadata_id, diverse_cohorts$name)
meta <- read.csv(meta_file$path)

if (verbose) {
  cat("\n", diverse_cohorts$name, "\n")
  print(colnames(meta))
  print_qc(meta, spec)
}

meta_new <- harmonize(diverse_cohorts$name, meta, spec) |>
  mutate(filename = meta_file$name)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", diverse_cohorts$name, "\n")
validate_values(meta_new, spec)

df_list[["Diverse_Cohorts"]] <- meta_new


# NPS-AD -----------------------------------------------------------------------

nps <- studies$nps_ad
meta_file <- synapse_download(nps$metadata_id, nps$name)
meta <- read.csv(meta_file$path)

neuro_file <- synapse_download(nps$extra_metadata, nps$name)
neuropath <- read.csv(neuro_file$path) |>
  dplyr::rename(individualID = IndividualID)

if (verbose) {
  cat("\n", nps$name, "\n")
  print(colnames(meta))
  print_qc(meta, spec)

  print_qc(neuropath, spec,
           column_renames = list(
             amyCerad = "CERAD",
             Braak = "BRAAK_AD"
           )
  )
}

meta_new <- harmonize(nps$name, meta, spec, extra_metadata = neuropath) |>
  mutate(filename = meta_file$name)

if (verbose) {
  print_qc(meta_new, spec)
}

cat("\n", nps$name, "\n")
validate_values(meta_new, spec) # Cohort will fail until after harmonization

df_list[["NPS-AD"]] <- meta_new


# Merge all files into one data frame ------------------------------------------

meta_all <- deduplicate_studies(df_list, spec, verbose = FALSE)

meta_all <- meta_all |>
  fill_missing_ampad1.0_ids(spec) |> # For Diverse Cohorts
  updateADoutcome(spec) # For Diverse Cohorts

if (verbose) {
  cat("\nAll studies (de-duplicated)\n")
  print_qc(meta_all, spec)
}

cat("\nAll studies (de-duplicated)\n")
validate_values(meta_all, spec)

# Save de-duplicated data frame for step 2
saveRDS(meta_all, file.path("data", "tmp", "AMP1.0_DiverseCohorts_harmonized.rds"))


# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
new_files <- sapply(df_list, function(meta_old) {
  meta_new <- subset(meta_all, filename == unique(meta_old$filename)) |>
    select(all_of(colnames(meta_old))) |>
    # Remove the source file column we added
    select(-filename)

  new_file <- str_replace(
    basename(unique(meta_old$filename)),
    "_harmonized.csv",
    ".csv"
  )

  # Sorting for this data set gets changed during fill_missing_ampad1.0_ids() so
  # we sort it back.
  if (unique(meta_new$study) == "NPS-AD") {
    meta_new <- arrange(meta_new, individualID)
  }

  # Now remove the study variable since it's no longer needed
  meta_new <- select(meta_new, -study)

  return(write_metadata(meta_new, new_file))
})

# Upload to ADKP metadata space

upload_id <- config::get("upload_synID")

for (filename in new_files) {
  synapse_upload(filename, upload_id)
}
