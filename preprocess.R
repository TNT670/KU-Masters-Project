#!/usr/bin/env Rscript
# preprocess.R
# Parameterized preprocessing pipeline for cancer tissue microarray data.
#
# Reads raw Affymetrix .CEL files from GSE tarballs, applies frozen RMA
# normalization, maps probes to Entrez gene IDs, collapses duplicate probes,
# applies quantile normalization, and removes batch effects via ComBat.
# Output is a space-delimited samples x genes expression matrix.
#
# Usage:
#   # All-inclusive run (all GSEs combined into one training matrix):
#   Rscript preprocess.R --cancer breast_cancer
#
#   # Leave-one-out run (hold out one GSE as a test set):
#   Rscript preprocess.R --cancer breast_cancer --loo 42568
#
#   # Custom config path:
#   Rscript preprocess.R --cancer lung_cancer --config /path/to/config.yaml
#
# Output files (written to the cancer type's directory):
#   All-inclusive:  {output_prefix}.dat
#   LOO train set:  {output_prefix}_loo{gse}_train.dat
#   LOO test set:   {output_prefix}_loo{gse}_test.dat
#
# NOTE: Extracted tarball contents are cached under {cancer_dir}/.extracted/
#   and are not automatically deleted between runs. Remove this directory
#   manually to force re-extraction.


# == 0. Package availability check ============================================
# optparse and yaml are CRAN packages; install them automatically if missing.
# Bioconductor packages (affy, frma, biomaRt, WGCNA, preprocessCore, sva)
# must be pre-installed via BiocManager::install().

for (pkg in c("optparse", "yaml")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(affy)
  library(frma)
  library(biomaRt)
  library(WGCNA)
  library(preprocessCore)
  library(sva)
})


# == 1. CLI argument parsing ==================================================

opt_list <- list(
  make_option("--cancer", type="character", default=NULL,
              help="Cancer type key from config.yaml (e.g. breast_cancer) [required]"),
  make_option("--loo",    type="character", default=NULL,
              help="GSE number to hold out as test set (e.g. 42568). Omit for all-inclusive run."),
  make_option("--config", type="character", default="config.yaml",
              help="Path to YAML config file [default: config.yaml]"),
  make_option("--local-annot", action="store_true", default=FALSE,
              help=paste("Use the local hgu133plus2.db Bioconductor package for",
                         "probe-to-Entrez mapping instead of querying Ensembl via biomaRt.",
                         "Faster, reproducible, and works offline.",
                         "Requires: BiocManager::install('hgu133plus2.db').",
                         "Mappings may differ slightly from biomaRt due to different curation."))
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$cancer)) {
  stop("--cancer is required. Example: Rscript preprocess.R --cancer breast_cancer")
}


# == 2. Config loading and validation =========================================

cat("Loading config from:", opt$config, "\n")
cfg_all <- yaml.load_file(opt$config)

if (!(opt$cancer %in% names(cfg_all))) {
  stop(sprintf("Cancer type '%s' not found in config. Available: %s",
               opt$cancer, paste(names(cfg_all), collapse=", ")))
}

cfg           <- cfg_all[[opt$cancer]]
cancer_dir    <- cfg$dir
output_prefix <- file.path(cancer_dir, cfg$output_prefix)
all_datasets  <- cfg$datasets
loo_gse       <- opt$loo

if (!is.null(loo_gse)) {
  all_gse_ids <- sapply(all_datasets, `[[`, "gse")
  if (!(loo_gse %in% all_gse_ids)) {
    stop(sprintf("LOO GSE '%s' not found in config for '%s'.", loo_gse, opt$cancer))
  }
}


# == 3. Helper functions ======================================================

# Expand a GSM range string like "GSM1615577-GSM1615614" into the full vector
# of all GSM IDs in that range (inclusive on both ends).
expand_gsm_range <- function(range_str) {
  parts     <- strsplit(range_str, "-")[[1]]
  start_num <- as.integer(sub("GSM", "", parts[1]))
  end_num   <- as.integer(sub("GSM", "", parts[2]))
  paste0("GSM", seq(start_num, end_num))
}

# Given a config `exclude` list (individual IDs and/or range strings),
# return a character vector of all excluded GSM IDs.
expand_exclusions <- function(excl_list) {
  if (is.null(excl_list) || length(excl_list) == 0) return(character(0))
  unlist(lapply(excl_list, function(e) {
    # Range format: "GSM1234567-GSM1234580"
    if (grepl("^GSM[0-9]+-GSM[0-9]+$", e)) expand_gsm_range(e) else e
  }))
}

# Extract a tarball into a subdirectory of base_dir and return the path to the
# directory containing the .CEL files.
# Handles tarballs that place all files in a single nested subdirectory.
# Skips extraction if CEL files are already present.
extract_tar <- function(tar_path, base_dir) {
  gse_tag <- sub("_RAW\\.tar$", "", basename(tar_path))
  out_dir <- file.path(base_dir, gse_tag)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  # Resolve the directory containing the CEL files, handling the case
  # where the tarball places everything inside a single nested subdirectory.
  resolve_cel_dir <- function(d) {
    entries <- list.files(d, full.names=TRUE)
    subdirs <- entries[file.info(entries)$isdir]
    if (length(subdirs) == 1 &&
        length(list.files(subdirs[[1]], pattern="\\.CEL(\\.gz)?$",
                          ignore.case=TRUE)) > 0) {
      return(subdirs[[1]])
    }
    d
  }

  # If CEL files already exist from a prior run, skip extraction.
  if (length(list.files(out_dir, pattern="\\.CEL(\\.gz)?$",
                         recursive=TRUE, ignore.case=TRUE)) > 0) {
    cat("  Using cached extraction for", basename(tar_path), "\n")
    return(resolve_cel_dir(out_dir))
  }

  cat("  Extracting", basename(tar_path), "...\n")
  untar(tar_path, exdir=out_dir)
  resolve_cel_dir(out_dir)
}

# List .CEL / .CEL.gz files in dir_path, excluding any whose leading GSM ID
# appears in excluded_ids. Returns a character vector of full file paths.
#
# Non-CEL files (e.g. .txt, .xml, .chp) are silently dropped here rather
# than being passed to ReadAffy, which would error on them.
get_cel_files <- function(dir_path, excluded_ids) {
  all_files <- list.files(dir_path, full.names=FALSE)

  # Filter to .CEL or .CEL.gz only
  cel_files <- all_files[grepl("\\.CEL(\\.gz)?$", all_files, ignore.case=TRUE)]
  if (length(cel_files) == 0) {
    stop(paste("No .CEL or .CEL.gz files found in", dir_path))
  }

  # Extract the GSM sample ID from each filename.
  # Filenames vary across GSEs; some use underscore ("GSM125119_array.CEL.gz"),
  # others use period ("GSM124994.CEL.gz"). Matching on the leading GSM+digits
  # token handles both formats reliably.
  gsm_prefix  <- regmatches(cel_files, regexpr("^GSM[0-9]+", cel_files))
  keep        <- !(gsm_prefix %in% excluded_ids)

  n_excl <- sum(!keep)
  if (n_excl > 0) {
    cat(sprintf("  Excluded %d file(s) matching exclusion list\n", n_excl))
  }

  file.path(dir_path, cel_files[keep])
}

# Read a set of .CEL files, apply frozen RMA (fRMA) normalization, and return
# the result as a samples x probes matrix.
#
# fRMA normalizes each array against a pre-computed frozen reference rather
# than within the current batch. This means arrays from different studies can
# be combined directly without the normalization introducing batch artifacts.
frma_normalize <- function(cel_files) {
  raw  <- ReadAffy(filenames=cel_files)
  eset <- frma(raw)
  # exprs() returns a probes x samples matrix; transpose to samples x probes
  t(as.data.frame(exprs(eset)))
}


# == 4. Split datasets into training and LOO sets =============================

if (!is.null(loo_gse)) {
  train_ds <- Filter(function(d) d$gse != loo_gse, all_datasets)
  loo_ds   <- Filter(function(d) d$gse == loo_gse, all_datasets)[[1]]
  cat(sprintf("\nLOO mode: %d training datasets, holding out GSE%s as test set\n\n",
              length(train_ds), loo_gse))
} else {
  train_ds <- all_datasets
  loo_ds   <- NULL
  cat(sprintf("\nAll-inclusive mode: %d datasets\n\n", length(train_ds)))
}

# Tarballs are extracted here and contents stored in the given directory
extract_base <- file.path(cancer_dir, ".extracted")
dir.create(extract_base, showWarnings=FALSE, recursive=TRUE)


# == 5. fRMA normalization, training datasets =================================
cat("=== Phase 1: fRMA normalization (training) ===\n")

dat         <- vector("list", length(train_ds))
batch_sizes <- integer(length(train_ds))  # post-exclusion sample counts per GSE

for (i in seq_along(train_ds)) {
  ds <- train_ds[[i]]
  cat(sprintf("[%d/%d] GSE%s\n", i, length(train_ds), ds$gse))

  tar_path       <- file.path(cancer_dir, sprintf("GSE%s_RAW.tar", ds$gse))
  excl_ids       <- expand_exclusions(ds$exclude)
  cel_dir        <- extract_tar(tar_path, extract_base)
  cel_files      <- get_cel_files(cel_dir, excl_ids)
  batch_sizes[i] <- length(cel_files)
  cat(sprintf("  Samples after exclusions: %d\n", batch_sizes[i]))

  dat[[i]] <- frma_normalize(cel_files)
}


# == 6. Combine all training samples ==========================================

cat("\n=== Phase 2: Combining training samples ===\n")
# Stack per-GSE sample matrices row-wise into one combined samples x probes matrix.
eset   <- do.call(rbind, dat)
cnames <- colnames(eset)  # probe IDs
cat(sprintf("Combined matrix: %d samples x %d probes\n", nrow(eset), ncol(eset)))


# == 7. Probe ID to Entrez gene ID mapping ====================================

cat("\n=== Phase 3: Probe to Entrez gene mapping ===\n")
# The HG-U133 Plus 2.0 array uses Affymetrix probe IDs. These are mapped to
# Entrez gene IDs to get a biologically interpretable gene matrix.
#
# If the --local-annot flag is omitted, biomaRt queries the Ensembl database
# via the internet for the mappings.
#
# If --local-annot is used, the mappings are obtained through the
# hgu133plus2.db Bioconductor package (offline, faster, and reproducible).
# Mappings may differ slightly from Ensembl.

if (isTRUE(opt$local_annot)) {

  cat("Annotation method: hgu133plus2.db (local)\n")
  if (!requireNamespace("hgu133plus2.db", quietly=TRUE)) {
    stop("Package 'hgu133plus2.db' is not installed. Run: BiocManager::install('hgu133plus2.db')")
  }
  suppressPackageStartupMessages(library(hgu133plus2.db))

  # select() may return multiple Entrez IDs per probe (one-to-many mappings);
  # keep only the first match per probe to maintain one Entrez ID per probe row
  entrez_map <- AnnotationDbi::select(hgu133plus2.db,
                                      keys    = cnames,
                                      columns = "ENTREZID",
                                      keytype = "PROBEID")
  entrez_map <- entrez_map[!duplicated(entrez_map$PROBEID), ]
  dftmp <- data.frame(
    cnames        = cnames,
    entrezgene_id = entrez_map$ENTREZID[match(cnames, entrez_map$PROBEID)]
  )
  cat(sprintf("%d of %d probes have a local Entrez mapping\n",
              sum(!is.na(dftmp$entrezgene_id)), length(cnames)))

} else {

  cat("Annotation method: biomaRt (Ensembl)\n")
  ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  annot   <- getBM(
    attributes = c("affy_hg_u133_plus_2", "entrezgene_id"),
    filters    = "affy_hg_u133_plus_2",
    values     = cnames,
    mart       = ensembl
  )
  if (nrow(annot) == 0) {
    stop(paste(
      "biomaRt returned no probe-to-gene mappings.",
      "Check that the Ensembl connection succeeded and that",
      "'affy_hg_u133_plus_2' is still the correct filter name.",
      "Alternatively, use --local-annot to bypass Ensembl entirely."
    ))
  }
  cat(sprintf("biomaRt returned mappings for %d probe entries\n", nrow(annot)))

  idx_lookup <- match(cnames, annot$affy_hg_u133_plus_2)
  dftmp      <- data.frame(cnames,
                           annot[idx_lookup, c("affy_hg_u133_plus_2", "entrezgene_id")])

}

# Transpose to probes x samples, attach Entrez IDs as a column,
# then drop any probe that has no Entrez match (NA).
eset_t        <- as.data.frame(t(eset))
eset_t$entrez <- dftmp$entrezgene_id
eset_t        <- na.omit(eset_t)
entrez_ids    <- eset_t$entrez
eset_t        <- eset_t[, names(eset_t) != "entrez"]
cat(sprintf("%d probes remain after removing unmapped probes\n", nrow(eset_t)))


# == 8. Collapse duplicate probes to single genes =============================

cat("\n=== Phase 4: Collapsing probes to genes ===\n")
# Several probes may map to the same Entrez gene ID.
# collapseRows selects the single most representative probe per gene (the probe
# with the highest mean expression across all samples).
collapsed    <- collapseRows(datET    = eset_t,
                             rowID    = rownames(eset_t),
                             rowGroup = entrez_ids)
collapsed_df <- data.frame(collapsed$datETcollapsed)
gene_ids     <- rownames(collapsed_df)
cat(sprintf("%d unique Entrez genes\n", length(gene_ids)))


# == 9. Quantile normalization ================================================
cat("\n=== Phase 5: Quantile normalization ===\n")
# Forces all samples to the same empirical expression distribution.
# This removes intensity-level differences between arrays that aren't
# attributable to biology.
# normalize.quantiles drops row and column names; restore them manually.
df_norm           <- as.data.frame(normalize.quantiles(as.matrix(collapsed_df)))
colnames(df_norm) <- colnames(collapsed_df)  # sample IDs
rownames(df_norm) <- gene_ids                # Entrez gene IDs


# == 10. ComBat batch effect removal ==========================================
cat("\n=== Phase 6: ComBat batch correction ===\n")
# Each GSE in the training set constitutes one batch (batch IDs 1..N).
# Sample counts per batch are derived from tarball contents after exclusions
# are applied
batch_labels <- rep(seq_along(train_ds), times=batch_sizes)
cat("Batch sizes (post-exclusion):", paste(batch_sizes, collapse=", "), "\n")

# ComBat models and removes systematic technical variation between batches
# (differences in lab protocols, scanner calibration, processing date, etc.)
# while preserving biological signal.
combat_data           <- ComBat(dat=df_norm, batch=batch_labels)
combat_data           <- as.data.frame(combat_data)
rownames(combat_data) <- gene_ids  # restore rownames

# Transpose to samples x genes for output
combat_t <- as.data.frame(t(combat_data))


# == 11. Write training output ================================================
cat("\n=== Phase 7: Writing training output ===\n")

output_df       <- combat_t
output_df$Batch <- batch_labels  # appended for downstream validation

if (is.null(loo_gse)) {
  out_path <- paste0(output_prefix, ".dat")
} else {
  out_path <- sprintf("%s_loo%s_train.dat", output_prefix, loo_gse)
}

write.table(output_df, out_path, quote=FALSE)
cat(sprintf("Written:    %s\n", out_path))
cat(sprintf("Dimensions: %d samples x %d genes (+1 Batch column)\n",
            nrow(output_df), ncol(output_df) - 1))


# == 12. LOO test set processing (only when --loo is specified) ===============

if (!is.null(loo_gse)) {

  cat(sprintf("\n=== LOO Phase 1: fRMA normalization — GSE%s ===\n", loo_gse))

  loo_tar     <- file.path(cancer_dir, sprintf("GSE%s_RAW.tar", loo_gse))
  loo_excl    <- expand_exclusions(loo_ds$exclude)
  loo_cel_dir <- extract_tar(loo_tar, extract_base)
  loo_cel     <- get_cel_files(loo_cel_dir, loo_excl)
  loo_n       <- length(loo_cel)
  cat(sprintf("LOO samples after exclusions: %d\n", loo_n))

  loo_eset <- frma_normalize(loo_cel)


  cat("\n=== LOO Phase 2: Probe → Entrez mapping (using training annotation) ===\n")
  # Reuse dftmp from the training run; this ensures the test set is mapped to
  # the identical gene space as the training set, which is required for the
  # downstream ComBat harmonization and classifier inference.
  loo_eset_t        <- as.data.frame(t(loo_eset))
  loo_eset_t$entrez <- dftmp$entrezgene_id
  loo_eset_t        <- na.omit(loo_eset_t)
  loo_entrez        <- loo_eset_t$entrez
  loo_eset_t        <- loo_eset_t[, names(loo_eset_t) != "entrez"]

  # Check that the LOO set resolves to the same gene set as the training set
  # (a mismatch would indicate a platform inconsistency across datasets)
  if (!setequal(loo_entrez, entrez_ids)) {
    stop(paste(
      "Gene space mismatch between training set and LOO GSE ", loo_gse,
      "; ensure all datasets use the same Affymetrix array platform."
    ))
  }


  cat("\n=== LOO Phase 3: Collapse probes to genes ===\n")
  loo_collapsed    <- collapseRows(datET    = loo_eset_t,
                                   rowID    = rownames(loo_eset_t),
                                   rowGroup = loo_entrez)
  loo_collapsed_df <- data.frame(loo_collapsed$datETcollapsed)


  cat("\n=== LOO Phase 4: Quantile normalization ===\n")
  loo_norm           <- as.data.frame(normalize.quantiles(as.matrix(loo_collapsed_df)))
  colnames(loo_norm) <- colnames(loo_collapsed_df)
  # Restore gene IDs as row names and reorder to match the training gene order
  # so rbind aligns columns correctly later
  rownames(loo_norm) <- rownames(loo_collapsed_df)
  loo_norm           <- loo_norm[gene_ids, , drop=FALSE]


  cat("\n=== LOO Phase 5: ComBat harmonization against training reference ===\n")
  # Merge the training samples (labeled batch 0) with the LOO test samples.
  # ComBat with ref.batch=0 adjusts the test set to the training distribution
  # without modifying the training data; the reference batch is the anchor.
  # After correction, only the test rows are extracted.
  merged            <- rbind(combat_t, t(loo_norm))  # (train_n + loo_n) x genes

  loo_batch_num     <- length(train_ds) + 1
  harmonize_batches <- c(rep(0L, nrow(combat_t)), rep(loo_batch_num, loo_n))

  combat_combined   <- ComBat(dat       = t(merged),
                              batch     = harmonize_batches,
                              ref.batch = 0)
  combined_t <- as.data.frame(t(combat_combined))  # all_samples x genes

  # Extract test rows by position
  n_train <- nrow(combat_t)
  cdtest  <- combined_t[(n_train + 1):nrow(combined_t), ]


  cat("\n=== LOO Phase 6: Writing test output ===\n")
  loo_out <- sprintf("%s_loo%s_test.dat", output_prefix, loo_gse)
  write.table(cdtest, loo_out, quote=FALSE)
  cat(sprintf("Written:    %s\n", loo_out))
  cat(sprintf("Dimensions: %d samples x %d genes\n", nrow(cdtest), ncol(cdtest)))
}

cat("\nPreprocessing complete.\n")
