## Pseudobulk differential expression on the merged Xenium object (from merge_xenium_h5ads.R)
## Sums raw counts per animal x cell type, then runs DESeq2 separately for each cell type.
## Handles any number of genotypes: all genotypes present for a cell type are fit in ONE model
## (shared dispersion estimates), then each non-reference genotype is tested against the reference
## (e.g. Het vs Wt, Hom vs Wt), plus any extra pairwise contrasts you list (e.g. Hom vs Het).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---------------- User settings ----------------
rds_path        <- "FMR1_merged_xenium_sce.rds"
sample_col      <- "animalID"        # biological replicate (from the per-sample YAMLs)
celltype_col    <- "subclass_name"   # cell type annotation column
geno_col        <- "genotype"
geno_recode     <- c(Ht = "Het")     # fix label errors: c(wrong = "right", ...); c() for none
geno_drop       <- c("Het")          # genotypes to exclude entirely (applied after recode); c() for none
geno_merge      <- list(HomHemi = c("Hom", "Hemi"))   # combine levels: list(new = c(old1, old2)); list() for none
# X-linked: males Hemi + females Hom -> one complete-loss group
ref_level       <- "Wt"
test_levels     <- NULL              # NULL = every non-reference genotype found (e.g. Het and Hom); or c("Het", "Hom")
extra_contrasts <- list()            # additional comparisons, e.g. list(c("Hom", "Het")) for Hom vs Het
covars_extra    <- c("animalSex")    # adjusted for when it has >1 level in a cell type; c() for none
min_cells       <- 20                # drop pseudobulk profiles built from fewer cells
min_reps        <- 3                 # min animals per genotype for that genotype to be included
control_genes   <- NULL              # optional stable genes for size factors, e.g. c("Gapdh", "Actb")
out_dir         <- "DEG_results/FMR1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

## ---------------- Load ----------------
sce <- readRDS(rds_path)

# Fix mislabeled genotypes before anything else
if (length(geno_recode)) {
  g <- as.character(sce[[geno_col]])
  hit <- g %in% names(geno_recode)
  if (any(hit)) {
    fixed <- unique(as.data.frame(colData(sce), optional = TRUE)[hit, c(sample_col, geno_col)])
    message("Recoding genotype for: ", paste0(fixed[[sample_col]], " (", fixed[[geno_col]], ")", collapse = ", "))
    g[hit] <- unname(geno_recode[g[hit]])
    sce[[geno_col]] <- g
  }
}

# Merge genotype levels (e.g. Hom + Hemi -> HomHemi)
for (new_lv in names(geno_merge)) {
  g <- as.character(sce[[geno_col]])
  old_lv <- geno_merge[[new_lv]]
  absent <- setdiff(old_lv, g)
  if (length(absent)) message("Note: merge level(s) not found in data: ", paste(absent, collapse = ", "))
  message("Merging ", paste(intersect(old_lv, g), collapse = " + "), " -> ", new_lv)
  g[g %in% old_lv] <- new_lv
  sce[[geno_col]] <- g
}

# Drop genotypes entirely (e.g. Het)
if (length(geno_drop)) {
  drop_cells <- sce[[geno_col]] %in% geno_drop
  n_anim <- length(unique(sce[[sample_col]][drop_cells]))
  message("Dropping ", paste(geno_drop, collapse = ", "), ": ", n_anim, " animals, ", sum(drop_cells), " cells")
  sce <- sce[, !drop_cells]
}
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
if (is.null(test_levels)) test_levels <- setdiff(geno_levels, ref_level)
bad_lv <- setdiff(c(test_levels, unlist(extra_contrasts)), geno_levels)
if (length(bad_lv)) stop("Genotype level(s) not found: ", paste(bad_lv, collapse = ", "))
if (length(test_levels) == 0) stop("No non-reference genotypes found")
message("Testing vs ", ref_level, ": ", paste(test_levels, collapse = ", "),
        if (length(extra_contrasts)) paste0(" | extra: ", paste(vapply(extra_contrasts, paste, "", collapse = " vs "), collapse = ", ")) else "")

message("Animals per genotype:")
animals <- cd_all %>% distinct(across(all_of(c(sample_col, geno_col, covars_extra))))
print(animals %>% count(across(all_of(c(geno_col, covars_extra)))))

# Sex balance check: if a genotype has only one sex, its effect can't be separated from sex.
if (length(covars_extra)) {
  for (v in covars_extra) {
    per_geno <- tapply(animals[[v]], animals[[geno_col]], function(z) length(unique(z)))
    one <- names(per_geno)[per_geno < 2]
    if (length(one) && length(unique(animals[[v]])) > 1) {
      warning("Genotype(s) ", paste(one, collapse = ", "), " have only one level of ", v,
              ". Genotype and ", v, " effects are partly confounded; interpret with care.")
    }
  }
}
message("Cells per cell type:")
print(sort(table(cd_all[[celltype_col]], useNA = "ifany"), decreasing = TRUE))

## ---------------- Remove control features, unannotated cells, unused genotypes ----------------
keep <- !grepl("NegControl|Unassigned|Deprecated|BLANK|Intergenic", rownames(sce))
if ("feature_types" %in% colnames(rowData(sce))) {
  keep <- keep & rowData(sce)$feature_types %in% "Gene Expression"
}
message("Keeping ", sum(keep), " of ", nrow(sce), " features")
sce <- sce[keep, !is.na(sce[[celltype_col]]) & sce[[geno_col]] %in% c(ref_level, test_levels)]

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
  
  # genotypes with enough animals in this cell type
  n_all <- table(factor(cd_ct[[geno_col]], levels = c(ref_level, test_levels)))
  ok_lv <- names(n_all)[n_all >= min_reps]
  if (!ref_level %in% ok_lv || length(ok_lv) < 2) {
    message("Skipping ", ct, ": animals per genotype = ", paste(names(n_all), n_all, collapse = ", "))
    next
  }
  dropped <- setdiff(c(ref_level, test_levels), ok_lv)
  if (length(dropped)) message("  ", ct, ": excluding ", paste(dropped, collapse = ", "), " (< ", min_reps, " animals)")
  in_model <- cd_ct[[geno_col]] %in% ok_lv
  sub <- sub[, in_model]; cd_ct <- cd_ct[in_model, , drop = FALSE]
  
  # keep covariates only if they vary in this cell type (e.g. one sex only -> drop sex)
  use_cov <- covars_extra[vapply(covars_extra, function(v) length(unique(cd_ct[[v]])) > 1, logical(1))]
  design_ct <- reformulate(c(use_cov, geno_col))
  for (v in use_cov) cd_ct[[v]] <- droplevels(factor(cd_ct[[v]]))
  cd_ct[[geno_col]] <- relevel(factor(cd_ct[[geno_col]], levels = ok_lv), ref = ref_level)
  
  cts <- as.matrix(assay(sub, cnt_assay))
  cts <- cts[rowSums(cts >= 10) >= min_reps, , drop = FALSE]   # light filter; panels are small
  rownames(cd_ct) <- colnames(cts)
  
  dds <- DESeqDataSetFromMatrix(round(cts), cd_ct, design_ct)
  if (!is.null(control_genes)) {
    dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% control_genes)
  }
  dds <- DESeq(dds, quiet = TRUE)
  
  # every available genotype vs reference, plus extra pairs where both genotypes are in the model
  contrasts <- c(lapply(setdiff(ok_lv, ref_level), function(l) c(l, ref_level)),
                 Filter(function(p) all(p %in% ok_lv), extra_contrasts))
  
  for (cn in contrasts) {
    cname <- paste(cn[1], "vs", cn[2])
    res <- results(dds, contrast = c(geno_col, cn[1], cn[2]))
    df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(contrast = cname, cell_type = ct,
             n_test = sum(cd_ct[[geno_col]] == cn[1]), n_ref = sum(cd_ct[[geno_col]] == cn[2]),
             genotypes_in_model = paste(ok_lv, collapse = "/"),
             design = deparse(design_ct)) %>%
      arrange(padj)
    
    message(ct, " | ", cname, ": ", sum(df$padj < 0.05, na.rm = TRUE), " DEGs (padj < 0.05)")
    cdir <- file.path(out_dir, safe(cname)); dir.create(cdir, showWarnings = FALSE)
    write.csv(df, file.path(cdir, paste0("DE_", safe(ct), ".csv")), row.names = FALSE)
    res_all[[paste(ct, cname)]] <- df
  }
}

if (length(res_all) == 0) {
  stop("No cell type had >= ", min_reps, " animals in ", ref_level, " and another genotype with >= ",
       min_cells, " cells. Lower min_reps/min_cells or check the pre-flight tables above.")
}
res_all <- bind_rows(res_all)
write.csv(res_all, file.path(out_dir, "DE_all_celltypes_all_contrasts.csv"), row.names = FALSE)

## ---------------- Summary ----------------
summ <- res_all %>%
  group_by(contrast, cell_type, n_test, n_ref) %>%
  summarise(n_up   = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
            n_down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
            .groups = "drop")
write.csv(summ, file.path(out_dir, "DE_summary.csv"), row.names = FALSE)
print(summ, n = Inf)

## ---------------- Volcano per contrast (written to file) ----------------
for (cname in unique(res_all$contrast)) {
  p <- ggplot(res_all %>% filter(contrast == cname, !is.na(padj)),
              aes(log2FoldChange, -log10(pvalue), colour = padj < 0.05)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(title = cname, colour = "padj < 0.05") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_dir, paste0("volcano_", safe(cname), ".png")), p, width = 10, height = 8, dpi = 300)
}

## ---------------- Dose check: compare fold changes of two genotypes vs reference ----------------
# e.g. Het vs Wt (x) against Hom vs Wt (y): a gene-dosage effect puts points along a line steeper than 1.
vs_ref <- paste(test_levels, "vs", ref_level)
vs_ref <- vs_ref[vs_ref %in% res_all$contrast]
if (length(vs_ref) >= 2) {
  wide <- res_all %>%
    filter(contrast %in% vs_ref[1:2]) %>%
    select(gene, cell_type, contrast, log2FoldChange, padj) %>%
    pivot_wider(names_from = contrast, values_from = c(log2FoldChange, padj)) %>%
    filter(!is.na(.data[[paste0("log2FoldChange_", vs_ref[1])]]),
           !is.na(.data[[paste0("log2FoldChange_", vs_ref[2])]])) %>%
    mutate(sig = case_when(
      .data[[paste0("padj_", vs_ref[1])]] < 0.05 & .data[[paste0("padj_", vs_ref[2])]] < 0.05 ~ "both",
      .data[[paste0("padj_", vs_ref[1])]] < 0.05 ~ vs_ref[1],
      .data[[paste0("padj_", vs_ref[2])]] < 0.05 ~ vs_ref[2],
      TRUE ~ "neither"))
  write.csv(wide, file.path(out_dir, paste0("LFC_compare_", safe(vs_ref[1]), "__", safe(vs_ref[2]), ".csv")), row.names = FALSE)
  
  p2 <- ggplot(wide, aes(.data[[paste0("log2FoldChange_", vs_ref[1])]],
                         .data[[paste0("log2FoldChange_", vs_ref[2])]], colour = sig)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
    geom_point(size = 0.8, alpha = 0.7) +
    facet_wrap(~ cell_type) +
    labs(x = paste("log2FC", vs_ref[1]), y = paste("log2FC", vs_ref[2]), colour = "padj < 0.05") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_dir, paste0("LFC_compare_", safe(vs_ref[1]), "__", safe(vs_ref[2]), ".png")),
         p2, width = 10, height = 8, dpi = 300)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Results in ", out_dir)