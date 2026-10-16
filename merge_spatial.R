## Load per-sample Xenium h5ad files + per-sample metadata YAMLs and merge into one SingleCellExperiment
## Expects pairs in one folder:  <name>.h5ad  +  <name>.yaml (or .yml)
## (the layout collect_xenium_samples.sh produces)
## Output: merged_xenium_sce.rds (input for xenium_pseudobulk_DE.R)

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(yaml)
  library(rhdf5)
})

## ---------------- User settings ----------------
h5ad_dir    <- "SHANK3"                    # folder with <name>.h5ad + <name>.yaml pairs
meta_fields <- NULL                      # NULL = keep all YAML fields; or e.g. c("animal", "sex", "genotype")
out_rds     <- "SHANK3_merged_xenium_sce.rds"

## ---------------- Find files ----------------
files <- list.files(h5ad_dir, pattern = "\\.h5ad$", full.names = TRUE)
names(files) <- sub("\\.h5ad$", "", basename(files))
if (length(files) == 0) stop("No .h5ad files in ", h5ad_dir)
message("Found ", length(files), " samples: ", paste(names(files), collapse = ", "))

yaml_for <- function(s) {
  cand <- file.path(h5ad_dir, paste0(s, c(".yaml", ".yml")))
  cand <- cand[file.exists(cand)]
  if (length(cand) == 0) stop("No metadata YAML for sample ", s, " (expected ", s, ".yaml or .yml)")
  cand[1]
}
yaml_files <- vapply(names(files), yaml_for, character(1))

## ---------------- Read per-sample metadata ----------------
# Nested fields are flattened with "." (e.g. slide: {id: X} -> slide.id); lists collapsed with ";".
# Returns a named character vector (one entry per field) for one sample.
read_one_yaml <- function(path, s) {
  y <- tryCatch(read_yaml(path), error = function(e) stop("Could not parse ", path, ":\n", conditionMessage(e)))
  if (!is.list(y) || length(y) == 0 || is.null(names(y))) {
    stop("YAML for ", s, " (", path, ") did not parse to key: value fields. Got: ",
         paste(capture.output(str(y)), collapse = " "))
  }
  flat <- lapply(y, function(v) {
    if (is.null(v)) return(NA_character_)
    if (is.list(v)) { u <- unlist(v); return(if (!is.null(names(u)) && all(nzchar(names(u)))) as.list(u) else paste(u, collapse = ";")) }
    paste(v, collapse = ";")
  })
  flat <- unlist(flat)                                  # joins nested names with "."
  flat <- setNames(as.character(flat), make.unique(names(flat)))
  names(flat)[names(flat) == "sample"] <- "yaml_sample"
  c(sample = s, yaml_file = basename(path), flat)
}

meta_list <- lapply(names(yaml_files), function(s) read_one_yaml(yaml_files[[s]], s))
names(meta_list) <- names(yaml_files)
message("Fields per YAML: ", paste(names(meta_list), lengths(meta_list) - 2, sep = "=", collapse = ", "))

# Base-R table: union of fields across samples, NA where a sample lacks a field; all character
fields <- unique(unlist(lapply(meta_list, names)))
meta <- as.data.frame(
  do.call(rbind, lapply(meta_list, function(v) setNames(v[fields], fields))),
  stringsAsFactors = FALSE, optional = TRUE
)
rownames(meta) <- NULL
if (!is.null(meta_fields)) {
  missing_f <- setdiff(meta_fields, names(meta))
  if (length(missing_f)) stop("Fields not in YAMLs: ", paste(missing_f, collapse = ", "),
                              "\nAvailable: ", paste(names(meta), collapse = ", "))
  meta <- meta[, c("sample", meta_fields), drop = FALSE]
}
na_cols <- names(meta)[colSums(is.na(meta)) > 0]
if (length(na_cols)) message("Fields missing in some YAMLs (NA): ", paste(na_cols, collapse = ", "))
print(meta)

## ---------------- Raw counts from adata.raw.X ----------------
# zellkonverter's R reader skips .raw, so read it straight from the HDF5 file.
h5_has_raw <- function(f) "raw" %in% rhdf5::h5ls(f, recursive = FALSE)$name

read_raw_X <- function(f) {
  enc   <- rhdf5::h5readAttributes(f, "raw/X")
  shape <- as.integer(enc$shape)                                  # cells, genes (Python order)
  d <- as.numeric(rhdf5::h5read(f, "raw/X/data"))
  i <- as.integer(rhdf5::h5read(f, "raw/X/indices"))
  p <- as.integer(rhdf5::h5read(f, "raw/X/indptr"))
  m <- if (enc$`encoding-type` == "csr_matrix") {
    Matrix::sparseMatrix(i = i + 1, p = p, x = d, dims = rev(shape))              # genes x cells
  } else {
    Matrix::t(Matrix::sparseMatrix(i = i + 1, p = p, x = d, dims = shape))        # csc -> genes x cells
  }
  gidx <- rhdf5::h5readAttributes(f, "raw/var")$`_index`
  cidx <- rhdf5::h5readAttributes(f, "obs")$`_index`
  dimnames(m) <- list(as.character(rhdf5::h5read(f, paste0("raw/var/", gidx))),
                      as.character(rhdf5::h5read(f, paste0("obs/", cidx))))
  m
}

attach_raw_counts <- function(x, f, s) {
  if (!h5_has_raw(f)) {
    message("  ", s, ": no raw.X; the DE script will rebuild counts from X")
    return(x)
  }
  raw <- read_raw_X(f)
  miss_g <- setdiff(rownames(x), rownames(raw)); miss_c <- setdiff(colnames(x), colnames(raw))
  if (length(miss_g) || length(miss_c)) {
    message("  ", s, ": raw.X is missing ", length(miss_g), " genes / ", length(miss_c),
            " cells present in X; skipping it (DE script will rebuild counts from X)")
    return(x)
  }
  raw <- raw[rownames(x), colnames(x)]
  if (any(raw@x != round(raw@x))) {
    message("  ", s, ": raw.X is not whole numbers; skipping it (DE script will rebuild counts from X)")
    return(x)
  }
  assay(x, "counts") <- raw   # X (log-normalized) is kept under its own name so all samples share it
  message("  ", s, ": counts from raw.X (", length(raw@x), " nonzero values)")
  x
}

## ---------------- Read each h5ad ----------------
sce_list <- lapply(names(files), function(s) {
  message("Reading ", s)
  x <- readH5AD(files[[s]], reader = "R")   # drop reader = "R" if your zellkonverter is older
  x <- attach_raw_counts(x, files[[s]], s)  # before renaming cells, so IDs match raw.X
  x$sample <- s
  colnames(x) <- paste(s, colnames(x), sep = "_")   # Xenium cell IDs repeat across samples
  for (v in names(colData(x))) if (is.factor(x[[v]])) x[[v]] <- as.character(x[[v]])
  metadata(x) <- list()                              # drop uns; often conflicts and wastes memory
  x
})
names(sce_list) <- names(files)

## ---------------- Harmonize features, assays, colData ----------------
genes     <- Reduce(intersect, lapply(sce_list, rownames))
all_genes <- Reduce(union,     lapply(sce_list, rownames))
if (length(genes) < length(all_genes)) {
  warning(length(all_genes) - length(genes), " features not shared by all samples were dropped, e.g. ",
          paste(head(setdiff(all_genes, genes), 20), collapse = ", "))
}

assays_keep <- Reduce(intersect, lapply(sce_list, assayNames))
no_counts <- names(sce_list)[!vapply(sce_list, function(x) "counts" %in% assayNames(x), logical(1))]
if (length(no_counts) && length(no_counts) < length(sce_list)) {
  message("No raw.X counts for: ", paste(no_counts, collapse = ", "),
          " -> 'counts' dropped for all samples; the DE script will rebuild counts from X")
}
rd_keep     <- intersect(Reduce(intersect, lapply(sce_list, reducedDimNames)), "spatial")
cd_cols     <- Reduce(union, lapply(sce_list, function(x) names(colData(x))))
row_info    <- rowData(sce_list[[1]])[genes, , drop = FALSE]
message("Assays kept: ", paste(assays_keep, collapse = ", "),
        " | reducedDims kept: ", paste(rd_keep, collapse = ", "))

sce_list <- lapply(sce_list, function(x) {
  x <- x[genes, ]
  assays(x)      <- assays(x)[assays_keep]
  reducedDims(x) <- reducedDims(x)[rd_keep]   # per-sample PCA/UMAP aren't comparable; keep spatial coords
  rowData(x)     <- NULL
  for (m in setdiff(cd_cols, names(colData(x)))) x[[m]] <- NA
  colData(x) <- colData(x)[, cd_cols]
  x
})

## ---------------- Merge ----------------
sce <- do.call(cbind, unname(sce_list))
rowData(sce) <- row_info
rm(sce_list); invisible(gc())

## ---------------- Attach sample metadata ----------------
clash <- setdiff(intersect(names(colData(sce)), names(meta)), "sample")
if (length(clash)) {
  message("Replacing colData columns with values from YAML: ", paste(clash, collapse = ", "))
  for (v in clash) sce[[v]] <- NULL
}
cd <- left_join(as.data.frame(colData(sce), optional = TRUE), meta, by = "sample")
colData(sce) <- DataFrame(cd, row.names = colnames(sce), check.names = FALSE)

## ---------------- Check + save ----------------
print(dim(sce))
print(as.data.frame(colData(sce)) %>% count(sample, across(any_of(c("animal", "sex", "genotype")))))
saveRDS(sce, out_rds)
message("Saved ", out_rds)