## Compare Xenium (spatial) pseudobulk DE with previous snRNA-seq DE
## For each matched snRNA <-> spatial cell type pair, on genes tested in both:
##   - Spearman correlation of log2 fold changes (uses all shared genes; most informative when DEGs are few)
##   - direction agreement: of the snRNA DEGs, how many change the same way in spatial (sign test)
##   - DEG overlap (Fisher test; universe = shared genes)
## Outputs CSVs + plots (written to file).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

## ---------------- User settings ----------------
spatial_file     <- "spatial_DEGs/DEG_results/FMR1/DE_all_celltypes_all_contrasts.csv"  # from xenium_pseudobulk_DE.R
spatial_contrast <- "HomHemi vs Wt"            # e.g. "Hom or Hemi vs Wt" 
# Spatial results use Ensembl IDs; snRNA uses gene symbols. Conversion source, first one set wins:
#   1) spatial_gene_table: your ID <-> name file (.csv, .tsv, .txt, optionally .gz; header optional)
spatial_gene_table    <- "spatial_DEGs/BTFTBJ_mBrain_480g_gene_list.csv"
gene_table_id_col     <- "Ensembl ID"    # column with Ensembl IDs; NULL = auto-detect (column of ENSMUSG... values)
gene_table_symbol_col <- "Gene"    # column with gene names; NULL = first other text column

# snRNA DE results: EITHER a folder of per-cell-type CSVs (cell type taken from the file name)
# OR one combined CSV with a cell type column.
snrna_path          <- "sexbygeno/DEGs/FMR1"
snrna_file_pattern  <- "\\.csv$"
snrna_celltype_from_file <- function(f) sub("\\.csv$", "", basename(f))   # edit if files are e.g. "edgeR_<celltype>.csv"
snrna_celltype_col  <- "cell_type"         # only used for a combined CSV
snrna_cols <- c(gene = "gene", lfc = "logFC", p = "PValue", fdr = "FDR")  # edgeR topTags names
# if no gene column, the first column (row names) is used
snrna_flip_sign <- FALSE   # TRUE if snRNA logFC is WT relative to Het (WT - Het); spatial is Het - Wt.
# Check: Myt1l itself should be DOWN in Het (see the printed check below)

# Cell type matching
# NULL = automatic: names are matched after removing the snRNA contrast suffix and treating
#        spaces, "/", "-", "_" as equivalent ("003 L5/6 IT TPE-ENT Glut" <-> "003-L5-6-IT-TPE-ENT-Glut_WT_vs_Het")
# Or a CSV with columns 'snrna' and 'spatial' for manual pairs (e.g. different taxonomies).
celltype_map_file <- NULL
snrna_name_suffix <- "_[^_]+_vs_[^_]+$"    # regex removed from snRNA names before matching (e.g. "_WT_vs_Het")

fdr_cut <- 0.05
out_dir <- "DE_comparison/FMR1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------- Read spatial ----------------
sp_all <- read_csv(spatial_file, show_col_types = FALSE)
if (!"contrast" %in% names(sp_all)) {
  # Older single-comparison DE output: every row is the one contrast that was run
  message("No 'contrast' column in ", basename(spatial_file), "; treating all rows as '", spatial_contrast, "'")
  sp_all$contrast <- spatial_contrast
}
if (!spatial_contrast %in% sp_all$contrast) {
  stop("Contrast '", spatial_contrast, "' not found. Available: ", paste(unique(sp_all$contrast), collapse = ", "))
}
## ---------------- Ensembl IDs -> gene symbols ----------------

sp_all$gene <- as.character(sp_all$gene)
is_ens <- mean(grepl("^ENSMUSG", sp_all$gene)) > 0.5
if (is_ens) {
  ens <- sub("\\.\\d+$", "", sp_all$gene)    # drop version suffix if present
  if (!is.null(spatial_gene_table)) {
    if (!file.exists(spatial_gene_table)) stop("spatial_gene_table not found: ", spatial_gene_table)
    sep <- if (grepl("\\.csv(\\.gz)?$", spatial_gene_table)) "," else "\t"
    tb <- read.delim(spatial_gene_table, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE, quote = "\"")
    if (any(grepl("^ENSMUSG", names(tb)))) {        # no header row (e.g. 10x features.tsv)
      tb <- read.delim(spatial_gene_table, sep = sep, header = FALSE, stringsAsFactors = FALSE, quote = "\"")
    }
    id_col <- if (!is.null(gene_table_id_col)) gene_table_id_col else
      names(tb)[which.max(vapply(tb, function(v) mean(grepl("^ENSMUSG", v)), numeric(1)))]
    sym_col <- if (!is.null(gene_table_symbol_col)) gene_table_symbol_col else
      setdiff(names(tb)[vapply(tb, is.character, logical(1))], id_col)[1]
    if (!all(c(id_col, sym_col) %in% names(tb)) || is.na(sym_col)) {
      stop("Could not find ID/name columns in ", spatial_gene_table, ". Columns: ", paste(names(tb), collapse = ", "),
           "\nSet gene_table_id_col and gene_table_symbol_col.")
    }
    message("Gene table: IDs from '", id_col, "', names from '", sym_col, "'. First rows:")
    print(head(tb[, c(id_col, sym_col)], 3))
    id2sym <- setNames(as.character(tb[[sym_col]]), sub("\\.\\d+$", "", as.character(tb[[id_col]])))
  } else if (!is.null(spatial_gene_h5ad)) {
    if (!file.exists(spatial_gene_h5ad)) stop("spatial_gene_h5ad not found: ", spatial_gene_h5ad)
    idx_col <- rhdf5::h5readAttributes(spatial_gene_h5ad, "var")$`_index`
    v <- list(index = read_h5_col(spatial_gene_h5ad, paste0("var/", idx_col)))
    for (cl in c("gene", "ensembl")) {
      v[[cl]] <- tryCatch(read_h5_col(spatial_gene_h5ad, paste0("var/", cl)), error = function(e) NULL)
    }
    if (is.null(v$gene)) stop("No 'gene' column in var of ", spatial_gene_h5ad, "; set spatial_gene_h5ad <- NULL to use org.Mm.eg.db")
    ens_src <- if (mean(grepl("^ENSMUSG", v$index)) > 0.5) v$index else v$ensembl
    if (is.null(ens_src)) stop("Could not find Ensembl IDs in var of ", spatial_gene_h5ad)
    id2sym <- setNames(v$gene, sub("\\.\\d+$", "", ens_src))
    message("Converting Ensembl IDs to symbols using the panel table in ", basename(spatial_gene_h5ad))
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Install org.Mm.eg.db or set spatial_gene_h5ad")
    id2sym <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = unique(ens),
                                    column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    message("Converting Ensembl IDs to symbols with org.Mm.eg.db")
  }
  sym <- unname(id2sym[ens])
  unmapped <- unique(sp_all$gene[is.na(sym) | sym == ""])
  if (length(unmapped)) message(length(unmapped), " Ensembl IDs had no symbol (kept as IDs; won't match snRNA): ",
                                paste(head(unmapped, 10), collapse = ", "))
  sp_all <- sp_all %>% mutate(ensembl = gene, gene = ifelse(is.na(sym) | sym == "", gene, sym))
  dup <- sp_all %>% count(cell_type, contrast, gene) %>% filter(n > 1)
  if (nrow(dup)) warning("Several Ensembl IDs map to the same symbol: ", paste(unique(dup$gene), collapse = ", "),
                         " -- these genes will appear more than once")
} else {
  sp_all$ensembl <- NA_character_
}

sp <- sp_all %>%
  filter(contrast == spatial_contrast) %>%
  transmute(gene = as.character(gene), ensembl, spatial = cell_type,
            lfc_sp = log2FoldChange, p_sp = pvalue, fdr_sp = padj)

## ---------------- Read snRNA ----------------
std_cols <- function(d, src) {
  gcol <- if (snrna_cols[["gene"]] %in% names(d)) snrna_cols[["gene"]] else names(d)[1]
  miss <- setdiff(snrna_cols[c("lfc", "p", "fdr")], names(d))
  if (length(miss)) stop("In ", src, ", column(s) not found: ", paste(miss, collapse = ", "),
                         "\nAvailable: ", paste(names(d), collapse = ", "), "\nEdit snrna_cols.")
  tibble(gene   = as.character(d[[gcol]]),
         lfc_sn = as.numeric(d[[snrna_cols[["lfc"]]]]),
         p_sn   = as.numeric(d[[snrna_cols[["p"]]]]),
         fdr_sn = as.numeric(d[[snrna_cols[["fdr"]]]]))
}

if (dir.exists(snrna_path)) {
  fs <- list.files(snrna_path, pattern = snrna_file_pattern, full.names = TRUE)
  if (!length(fs)) stop("No files matching ", snrna_file_pattern, " in ", snrna_path)
  sn <- bind_rows(lapply(fs, function(f) {
    std_cols(read_csv(f, show_col_types = FALSE), f) %>% mutate(snrna = snrna_celltype_from_file(f))
  }))
} else {
  d <- read_csv(snrna_path, show_col_types = FALSE)
  if (!snrna_celltype_col %in% names(d)) stop("No column '", snrna_celltype_col, "' in ", snrna_path)
  sn <- std_cols(d, snrna_path) %>% mutate(snrna = as.character(d[[snrna_celltype_col]]))
}

if (snrna_flip_sign) {
  sn$lfc_sn <- -sn$lfc_sn
  message("snRNA log2FC sign flipped so positive = up in mutant")
}
sanity_gene <- "Myt1l"
if (sanity_gene %in% sn$gene) {
  message(sanity_gene, " snRNA log2FC by cell type (should be negative = down in Het):")
  print(sn %>% filter(gene == sanity_gene) %>% select(snrna, lfc_sn, fdr_sn), n = Inf)
}

# Always write both cell type lists -- handy for building celltype_map.csv
write_csv(tibble(snrna = sort(unique(sn$snrna))), file.path(out_dir, "celltypes_snrna.csv"))
write_csv(sp %>% count(spatial, name = "n_genes_tested") %>% arrange(spatial), file.path(out_dir, "celltypes_spatial.csv"))

## ---------------- Panel coverage: how many snRNA DEGs could the Xenium panel even detect? ----------------
panel_genes <- unique(sp_all$gene)
coverage <- sn %>%
  filter(!is.na(fdr_sn), fdr_sn < fdr_cut) %>%
  group_by(snrna) %>%
  summarise(n_snrna_DEGs = n(), n_on_panel = sum(gene %in% panel_genes),
            DEGs_on_panel = paste(sort(gene[gene %in% panel_genes]), collapse = ";"), .groups = "drop")
write_csv(coverage, file.path(out_dir, "snrna_DEGs_on_xenium_panel.csv"))
message("snRNA DEGs on the Xenium panel: ", sum(coverage$n_on_panel), " of ", sum(coverage$n_snrna_DEGs))
if (sum(sn$gene %in% panel_genes) == 0) {
  stop("No snRNA genes match spatial genes -- check gene ID types (symbols vs Ensembl) and capitalization.")
}

## ---------------- Cell type map ----------------
if (!is.null(celltype_map_file)) {
  ct_map <- read_csv(celltype_map_file, show_col_types = FALSE)
  if (!all(c("snrna", "spatial") %in% names(ct_map))) stop(celltype_map_file, " needs columns 'snrna' and 'spatial'")
  ct_map <- ct_map %>% select(snrna, spatial) %>% distinct()
  bad_sn <- setdiff(ct_map$snrna, sn$snrna); bad_sp <- setdiff(ct_map$spatial, sp$spatial)
  if (length(bad_sn)) message("Map snRNA names with no results: ", paste(bad_sn, collapse = ", "))
  if (length(bad_sp)) message("Map spatial names with no results (not tested?): ", paste(bad_sp, collapse = ", "))
} else {
  norm_name <- function(x) {
    k <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
    gsub("^-+|-+$", "", k)
  }
  sn_names <- tibble(snrna   = sort(unique(sn$snrna))) %>% mutate(key = norm_name(sub(snrna_name_suffix, "", snrna)))
  sp_names <- tibble(spatial = sort(unique(sp$spatial))) %>% mutate(key = norm_name(spatial))
  dup <- c(sn_names$key[duplicated(sn_names$key)], sp_names$key[duplicated(sp_names$key)])
  if (length(dup)) warning("Names that collapse to the same key (check these pairs): ", paste(unique(dup), collapse = ", "))
  
  ct_map <- inner_join(sn_names, sp_names, by = "key") %>% select(snrna, spatial)
  unmatched_sp <- setdiff(sp_names$spatial, ct_map$spatial)
  unmatched_sn <- setdiff(sn_names$snrna, ct_map$snrna)
  message("Auto-matched ", nrow(ct_map), " cell types")
  if (length(unmatched_sp)) message("  Spatial with no snRNA match: ", paste(unmatched_sp, collapse = "; "))
  if (length(unmatched_sn)) message("  snRNA with no spatial match: ", paste(unmatched_sn, collapse = "; "))
}
write_csv(ct_map, file.path(out_dir, "celltype_map_used.csv"))
if (nrow(ct_map) == 0) {
  stop("No matching cell types. See celltypes_snrna.csv and celltypes_spatial.csv in ", out_dir,
       " and make a celltype_map.csv (columns: snrna, spatial).")
}

## ---------------- Join on shared genes ----------------
joined <- ct_map %>%
  inner_join(sn, by = "snrna") %>%
  inner_join(sp, by = c("spatial", "gene")) %>%
  mutate(sig_sn = !is.na(fdr_sn) & fdr_sn < fdr_cut,
         sig_sp = !is.na(fdr_sp) & fdr_sp < fdr_cut,
         category = case_when(sig_sn & sig_sp ~ "both",
                              sig_sn ~ "snRNA only",
                              sig_sp ~ "spatial only",
                              TRUE ~ "neither"),
         same_direction = sign(lfc_sn) == sign(lfc_sp),
         pair = paste(snrna, "|", spatial))
write_csv(joined, file.path(out_dir, "genes_compared.csv"))

## ---------------- Per-pair statistics ----------------
pair_stats <- function(d, key) {
  n <- nrow(d)
  rho <- rho_p <- NA_real_
  if (n >= 5) {
    ct <- suppressWarnings(cor.test(d$lfc_sn, d$lfc_sp, method = "spearman", exact = FALSE))
    rho <- unname(ct$estimate); rho_p <- ct$p.value
  }
  deg <- d[d$sig_sn & !is.na(d$lfc_sp), ]
  n_sn <- nrow(deg); same <- sum(deg$same_direction)
  both <- sum(d$sig_sn & d$sig_sp); sn_only <- sum(d$sig_sn & !d$sig_sp)
  sp_only <- sum(!d$sig_sn & d$sig_sp); neither <- sum(!d$sig_sn & !d$sig_sp)
  tibble(
    n_genes_shared        = n,
    spearman_rho          = rho,
    rho_p                 = rho_p,
    n_deg_snrna           = sum(d$sig_sn),
    n_deg_spatial         = sum(d$sig_sp),
    n_deg_both            = both,
    snrna_degs_same_dir   = same,
    frac_same_dir         = if (n_sn) same / n_sn else NA_real_,
    sign_test_p           = if (n_sn) binom.test(same, n_sn, 0.5, alternative = "greater")$p.value else NA_real_,
    overlap_fisher_p      = if (both + sn_only > 0 && both + sp_only > 0)
      fisher.test(matrix(c(both, sn_only, sp_only, neither), 2), alternative = "greater")$p.value
    else NA_real_,
    genes_deg_both        = paste(sort(d$gene[d$sig_sn & d$sig_sp]), collapse = ";")
  )
}

summ <- joined %>%
  group_by(snrna, spatial, pair) %>%
  group_modify(pair_stats) %>%
  ungroup() %>%
  mutate(rho_fdr = p.adjust(rho_p, "BH")) %>%
  arrange(desc(spearman_rho))
write_csv(summ, file.path(out_dir, "comparison_summary.csv"))
print(summ %>% select(pair, n_genes_shared, spearman_rho, rho_fdr, n_deg_snrna, frac_same_dir, n_deg_both), n = Inf)

## ---------------- Plots ----------------
cat_cols <- c(both = "firebrick", `snRNA only` = "steelblue", `spatial only` = "darkorange", neither = "grey75")

ann <- summ %>% transmute(pair, label = sprintf("rho = %.2f\nn = %d", spearman_rho, n_genes_shared))
n_pairs <- length(unique(joined$pair)); ncol_f <- ceiling(sqrt(n_pairs))

p1 <- ggplot(joined %>% arrange(factor(category, levels = rev(names(cat_cols)))),
             aes(lfc_sn, lfc_sp, colour = category)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 0.9, alpha = 0.8) +
  geom_text(data = joined %>% filter(category == "both"), aes(label = gene),
            size = 2.3, vjust = -0.6, show.legend = FALSE, check_overlap = TRUE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.15, size = 2.6) +
  scale_colour_manual(values = cat_cols) +
  facet_wrap(~ pair, scales = "free", ncol = ncol_f) +
  labs(x = "snRNA-seq log2FC", y = paste("Xenium log2FC,", spatial_contrast),
       colour = paste0("FDR < ", fdr_cut)) +
  theme_bw(base_size = 9)
ggsave(file.path(out_dir, "LFC_scatter_by_celltype.png"), p1,
       width = max(6, 3.2 * ncol_f), height = max(5, 3 * ceiling(n_pairs / ncol_f)), dpi = 300, limitsize = FALSE)

p2 <- ggplot(summ %>% filter(!is.na(spearman_rho)),
             aes(reorder(pair, spearman_rho), spearman_rho, colour = rho_fdr < 0.05, size = n_genes_shared)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point() +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick")) +
  coord_flip() +
  labs(x = NULL, y = "Spearman rho of log2FC (snRNA vs Xenium)", colour = "BH FDR < 0.05", size = "genes") +
  theme_bw(base_size = 9)
ggsave(file.path(out_dir, "LFC_correlation_summary.png"), p2,
       width = 7, height = max(3, 0.3 * nrow(summ) + 1.5), dpi = 300, limitsize = FALSE)

message("Results in ", out_dir)
