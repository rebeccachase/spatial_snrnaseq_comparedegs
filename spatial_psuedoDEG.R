## Pseudobulk differential expression on the merged Xenium object (from merge_xenium_h5ads.R)
## Sums raw counts per animal x cell type, then runs DESeq2 separately for each cell type.
#this is for het v wt only 
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

## ---------------- User settings ----------------
rds_path      <- "NF1_merged_xenium_sce.rds"
sample_col    <- "animalID"       # biological replicate (from the per-sample YAMLs)
celltype_col  <- "subclass_name"      # cell type annotation column in the h5ad obs -- check the pre-flight output
geno_col      <- "genotype"
ref_level     <- "Wt"
test_level    <- NULL             # NULL = use the only non-reference genotype; or set e.g. "Het"
covars_extra  <- c("animalSex")   # adjusted for when it has >1 level in a cell type; c() for none
min_cells     <- 20               # drop pseudobulk profiles built from fewer cells
min_reps      <- 3                # min animals per genotype to test a cell type
control_genes <- NULL             # optional stable genes for size factors, e.g. c("Gapdh", "Actb")
out_dir       <- "DEG_results/NF1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------- Load ----------------
sce <- readRDS(rds_path)

#have to do this for NF1 because of errors in genotype
#g <- as.character(sce$genotype)
#sce$genotype <- ifelse(g %in% "Ht", "Het", g)
#table(sce$genotype)   # check: no "Ht" left

cnt_assay <- if ("counts" %in% assayNames(sce)) "counts" else "X"
message("Using assay: ", cnt_assay, " | ", nrow(sce), " genes x ", ncol(sce), " cells")

chk <- as.matrix(assay(sce, cnt_assay)[, seq_len(min(500, ncol(sce)))])
if (any(chk != round(chk))) stop("Assay '", cnt_assay, "' is not raw counts. Point to the raw layer.")

## ---------------- Pre-flight checks ----------------
cd_all <- as.data.frame(colData(sce), optional = TRUE)
needed <- c(sample_col, celltype_col, geno_col, covars_extra)
missing_cols <- setdiff(needed, names(cd_all))
if (length(missing_cols)) {
  cand <- names(cd_all)[vapply(cd_all, function(v) (is.character(v) || is.factor(v)) &&
                                 length(unique(v)) < 200, logical(1))]
  stop("Missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nCategorical columns available: ", paste(cand, collapse = ", "))
}

geno_levels <- sort(unique(as.character(cd_all[[geno_col]])))
message("Genotype levels: ", paste(geno_levels, collapse = ", "))
if (!ref_level %in% geno_levels) stop("Reference level '", ref_level, "' not in ", geno_col)
if (is.null(test_level)) {
  others <- setdiff(geno_levels, ref_level)
  if (length(others) != 1) stop("Several non-reference genotypes (", paste(others, collapse = ", "),
                                "); set test_level")
  test_level <- others
}
message("Testing ", test_level, " vs ", ref_level)

message("Animals per genotype:")
print(cd_all %>% distinct(across(all_of(c(sample_col, geno_col, covars_extra)))) %>%
        count(across(all_of(c(geno_col, covars_extra)))))
message("Cells per cell type:")
print(sort(table(cd_all[[celltype_col]], useNA = "ifany"), decreasing = TRUE))

## ---------------- Remove control features and unannotated cells ----------------
keep <- !grepl("NegControl|Unassigned|Deprecated|BLANK|Intergenic", rownames(sce))
if ("feature_types" %in% colnames(rowData(sce))) {
  keep <- keep & rowData(sce)$feature_types %in% "Gene Expression"
}
message("Keeping ", sum(keep), " of ", nrow(sce), " features")
sce <- sce[keep, !is.na(sce[[celltype_col]])]

## ---------------- Sample-level metadata ----------------
meta_cols <- c(sample_col, geno_col, covars_extra)
sample_meta <- as.data.frame(colData(sce), optional = TRUE)[, meta_cols, drop = FALSE] %>% distinct()
if (anyDuplicated(sample_meta[[sample_col]])) {
  stop("Genotype/covariates vary within an animal ID; check ", paste(meta_cols[-1], collapse = ", "))
}

## ---------------- Pseudobulk ----------------
pb <- aggregateAcrossCells(
  sce,
  ids = colData(sce)[, c(sample_col, celltype_col)],
  use.assay.type = cnt_assay,
  statistics = "sum"
)
cd <- as.data.frame(colData(pb), optional = TRUE)[, c(sample_col, celltype_col, "ncells")]
cd <- left_join(cd, sample_meta, by = sample_col)   # preserves row order
colData(pb) <- DataFrame(cd, check.names = FALSE)
colnames(pb) <- paste(pb[[sample_col]], pb[[celltype_col]], sep = "__")
write.csv(cd, file.path(out_dir, "pseudobulk_sample_table.csv"), row.names = FALSE)

## ---------------- DE per cell type ----------------
res_all <- list()

for (ct in sort(unique(pb[[celltype_col]]))) {
  sub <- pb[, pb[[celltype_col]] == ct & pb$ncells >= min_cells]
  cd_ct <- as.data.frame(colData(sub), optional = TRUE)
  
  n_grp <- table(factor(cd_ct[[geno_col]], levels = c(test_level, ref_level)))
  if (any(n_grp < min_reps)) {
    message("Skipping ", ct, ": animals per genotype = ", paste(names(n_grp), n_grp, collapse = ", "))
    next
  }
  
  # keep covariates only if they vary in this cell type (e.g. one sex only -> drop sex)
  use_cov <- covars_extra[vapply(covars_extra, function(v) length(unique(cd_ct[[v]])) > 1, logical(1))]
  design_ct <- reformulate(c(use_cov, geno_col))
  for (v in c(use_cov, geno_col)) cd_ct[[v]] <- droplevels(factor(cd_ct[[v]]))
  cd_ct[[geno_col]] <- relevel(cd_ct[[geno_col]], ref = ref_level)
  
  cts <- as.matrix(assay(sub, cnt_assay))
  cts <- cts[rowSums(cts >= 10) >= min_reps, , drop = FALSE]   # light filter; panels are small
  rownames(cd_ct) <- colnames(cts)
  
  dds <- DESeqDataSetFromMatrix(round(cts), cd_ct, design_ct)
  if (!is.null(control_genes)) {
    dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% control_genes)
  }
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c(geno_col, test_level, ref_level))
  
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(cell_type = ct, n_animals = ncol(dds),
           n_test = n_grp[[test_level]], n_ref = n_grp[[ref_level]],
           design = deparse(design_ct)) %>%
    arrange(padj)
  
  message(ct, ": ", sum(df$padj < 0.05, na.rm = TRUE), " DEGs (padj < 0.05), design ", deparse(design_ct))
  write.csv(df, file.path(out_dir, paste0("DE_", gsub("[^A-Za-z0-9]+", "_", ct), ".csv")), row.names = FALSE)
  res_all[[ct]] <- df
}

if (length(res_all) == 0) {
  stop("No cell type had >= ", min_reps, " animals per genotype with >= ", min_cells,
       " cells. Lower min_reps/min_cells or check the pre-flight tables above.")
}
res_all <- bind_rows(res_all)
write.csv(res_all, file.path(out_dir, "DE_all_celltypes.csv"), row.names = FALSE)

## ---------------- Summary + volcano (written to file) ----------------
summ <- res_all %>%
  group_by(cell_type, n_test, n_ref) %>%
  summarise(n_up   = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
            n_down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
            .groups = "drop")
write.csv(summ, file.path(out_dir, "DE_summary.csv"), row.names = FALSE)
print(summ)

p <- ggplot(res_all %>% filter(!is.na(padj)),
            aes(log2FoldChange, -log10(pvalue), colour = padj < 0.05)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(title = paste(test_level, "vs", ref_level), colour = "padj < 0.05") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "volcano_by_celltype.png"), p, width = 10, height = 8, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Results in ", out_dir)
