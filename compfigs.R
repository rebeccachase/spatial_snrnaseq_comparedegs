## Figures summarizing the Xenium vs snRNA-seq DE comparison, for each mouse line + a cross-line overview
## Reads the outputs of compare_spatial_snrna_DE.R (run once per line, each with its own out_dir);
## writes PNGs + a multi-page PDF per line into <line folder>/figures, and overview figures into overview_dir.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

## ---------------- Settings ----------------
# One entry per comparison folder (name = label used in titles). Lines with Het AND Hom (e.g. SHANK3)
# need one compare run per contrast, each with its own out_dir.
comparison_dirs <- c(
  FMR1       = "DE_comparison/FMR1",
  DNMT3A     = "DE_comparison/DNMT3A",
  MYT1L      = "DE_comparison/MYT1L",
  SHANK3_Het = "DE_comparison/SHANK3/Het_v_WT",
  SHANK3_Hom = "DE_comparison/SHANK3/Hom_v_WT"
)
fdr_cut      <- 0.05
overview_dir <- "DE_comparison/overview"
dir.create(overview_dir, showWarnings = FALSE, recursive = TRUE)

cell_class <- function(x) sub(".*\\s", "", x)   # last word of Allen name: Glut, Gaba, NN, IMN, ...

## ---------------- Per-line figures ----------------
plot_line <- function(in_dir, line) {
  out_dir <- file.path(in_dir, "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ## ---------------- Data ----------------
  cell_class <- function(x) sub(".*\\s", "", x)   # last word of Allen name: Glut, Gaba, NN, IMN, ...
  
  summ <- read_csv(file.path(in_dir, "comparison_summary.csv"), show_col_types = FALSE) %>%
    filter(!is.na(spearman_rho)) %>%
    mutate(class = cell_class(spatial),
           sig   = ifelse(!is.na(rho_fdr) & rho_fdr < fdr_cut, paste0("FDR < ", fdr_cut), "n.s."))
  genes <- read_csv(file.path(in_dir, "genes_compared.csv"), show_col_types = FALSE) %>%
    mutate(class = cell_class(spatial))
  # results from the older compare script (no sex/group split) have no snrna_group column
  if (!"snrna_group" %in% names(summ))  summ$snrna_group  <- "all"
  if (!"snrna_group" %in% names(genes)) genes$snrna_group <- "all"
  groups <- sort(unique(summ$snrna_group))
  plots <- list()
  
  ## ---------------- A) Correlation per cell type ----------------
  ord <- summ %>% group_by(spatial) %>% summarise(m = mean(spearman_rho)) %>% arrange(m) %>% pull(spatial)
  plots$A_rho_by_celltype <- ggplot(summ, aes(spearman_rho, factor(spatial, levels = ord))) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(aes(colour = class, shape = sig, size = n_deg_snrna)) +
    scale_shape_manual(values = setNames(c(16, 1), c(paste0("FDR < ", fdr_cut), "n.s."))) +
    scale_size(range = c(1, 4)) +
    facet_wrap(~ snrna_group) +
    labs(x = "Spearman rho of log2FC (snRNA vs Xenium)", y = NULL, colour = "class",
         shape = "correlation", size = "snRNA DEGs\non panel",
         title = paste0(line, " - ", "Agreement per cell type")) +
    theme_bw(base_size = 8)
  attr(plots$A_rho_by_celltype, "dims") <- c(4 + 3 * length(groups), max(5, 0.11 * length(ord) + 2))
  
  ## ---------------- B) Distribution of correlations ----------------
  stats_B <- summ %>% group_by(snrna_group) %>%
    summarise(n = n(), pos = sum(spearman_rho > 0), med = median(spearman_rho),
              p = binom.test(pos, n, 0.5, alternative = "greater")$p.value, .groups = "drop") %>%
    mutate(label = sprintf("%d / %d cell types positive\nmedian rho = %.2f\nsign test p = %.2g", pos, n, med, p))
  plots$B_rho_distribution <- ggplot(summ, aes(spearman_rho)) +
    geom_histogram(bins = 30, fill = "grey70", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(data = stats_B, aes(xintercept = med), colour = "firebrick") +
    geom_text(data = stats_B, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2, size = 3) +
    facet_wrap(~ snrna_group) +
    labs(x = "Spearman rho (one value per cell type)", y = "cell types",
         title = paste0(line, " - ", "Are correlations shifted above zero?"), subtitle = "red = median; dashed = no agreement") +
    theme_bw(base_size = 10)
  attr(plots$B_rho_distribution, "dims") <- c(4 + 3 * length(groups), 4)
  
  ## ---------------- C) Direction agreement of snRNA DEGs ----------------
  degs <- genes %>% filter(sig_sn, !is.na(same_direction)) %>%
    mutate(direction = ifelse(same_direction, "same", "opposite"))
  if (nrow(degs)) {
    pooled <- degs %>% group_by(snrna_group) %>%
      summarise(n = n(), same = sum(same_direction),
                p = binom.test(same, n, 0.5, alternative = "greater")$p.value, .groups = "drop")
    plots$C_direction_agreement <- ggplot(degs, aes(class, fill = direction)) +
      geom_bar(position = "fill") +
      geom_text(stat = "count", aes(label = after_stat(count)), position = position_fill(vjust = 0.5), size = 2.6) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      scale_fill_manual(values = c(same = "steelblue", opposite = "grey75")) +
      facet_wrap(~ snrna_group) +
      labs(x = NULL, y = "fraction of snRNA DEGs", fill = "Xenium direction",
           title = paste0(line, " - ", "Do snRNA DEGs change the same way in Xenium?"),
           subtitle = paste(sprintf("%s: %d / %d same (%.0f%%), p = %.2g",
                                    pooled$snrna_group, pooled$same, pooled$n, 100 * pooled$same / pooled$n, pooled$p),
                            collapse = "    ")) +
      theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    attr(plots$C_direction_agreement, "dims") <- c(4 + 3 * length(groups), 4.5)
  }
  
  ## ---------------- D) Fold changes, all cell types pooled ----------------
  # each gene labeled once (at its largest Xenium change), with how many cell types it is DE in both
  lab_D <- genes %>% filter(category == "both") %>% group_by(snrna_group, gene) %>%
    summarise(n_ct = n(), i = which.max(abs(lfc_sp)), lfc_sn = lfc_sn[i], lfc_sp = lfc_sp[i], .groups = "drop") %>%
    mutate(label = ifelse(n_ct > 1, paste0(gene, " (", n_ct, " cell types)"), gene))
  lims_x <- quantile(genes$lfc_sn, c(0.002, 0.998), na.rm = TRUE)
  lims_y <- quantile(genes$lfc_sp, c(0.002, 0.998), na.rm = TRUE)
  plots$D_lfc_pooled <- ggplot(genes, aes(lfc_sn, lfc_sp)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(data = filter(genes, category == "neither"), colour = "grey85", size = 0.3, alpha = 0.4) +
    geom_point(data = filter(genes, category != "neither"), aes(colour = category), size = 0.9, alpha = 0.8) +
    geom_smooth(data = filter(genes, sig_sn), method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = 0.5) +
    geom_text(data = lab_D, aes(label = label), size = 2.3, vjust = -0.6, check_overlap = TRUE) +
    scale_colour_manual(values = c(both = "firebrick", `snRNA only` = "steelblue", `spatial only` = "darkorange")) +
    coord_cartesian(xlim = lims_x, ylim = lims_y) +
    facet_wrap(~ snrna_group) +
    labs(x = "snRNA-seq log2FC", y = "Xenium log2FC", colour = paste0("FDR < ", fdr_cut),
         title = paste0(line, " - ", "All gene x cell type fold changes"),
         subtitle = "grey = not DE in either; line = fit through snRNA DEGs") +
    theme_bw(base_size = 10)
  attr(plots$D_lfc_pooled, "dims") <- c(4 + 3.5 * length(groups), 4.5)
  
  ## ---------------- E) Group vs group (e.g. f-hom vs m-hemi) ----------------
  if (length(groups) >= 2) {
    wide <- summ %>% select(snrna_group, spatial, class, spearman_rho) %>%
      pivot_wider(names_from = snrna_group, values_from = spearman_rho) %>%
      filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]])) %>%
      mutate(diff = abs(.data[[groups[1]]] - .data[[groups[2]]]))
    plots$E_group_comparison <- ggplot(wide, aes(.data[[groups[1]]], .data[[groups[2]]], colour = class)) +
      geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
      geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
      geom_point(size = 1.8) +
      geom_text(data = slice_max(wide, diff, n = 6), aes(label = spatial), size = 2.4, vjust = -0.7,
                show.legend = FALSE, check_overlap = TRUE) +
      labs(x = paste("rho,", groups[1]), y = paste("rho,", groups[2]),
           title = paste0(line, " - agreement in ", groups[1], " vs ", groups[2]),
           subtitle = "one point per cell type; labeled = largest differences") +
      theme_bw(base_size = 10)
    attr(plots$E_group_comparison, "dims") <- c(6.5, 5.5)
  }
  
  ## ---------------- Write ----------------
  for (nm in names(plots)) {
    d <- attr(plots[[nm]], "dims")
    ggsave(file.path(out_dir, paste0(nm, ".png")), plots[[nm]], width = d[1], height = d[2], dpi = 300, limitsize = FALSE)
  }
  pdf(file.path(out_dir, "comparison_figures.pdf"), width = 11, height = 8.5)
  for (nm in names(plots)) print(plots[[nm]])
  dev.off()
  message(line, ": figures in ", out_dir)
  list(summ = summ %>% mutate(line = line), degs = degs %>% mutate(line = line))
}

found <- comparison_dirs[file.exists(file.path(comparison_dirs, "comparison_summary.csv"))]
missing <- setdiff(names(comparison_dirs), names(found))
if (length(missing)) message("Skipping (no comparison_summary.csv yet): ", paste(missing, collapse = ", "))
if (!length(found)) stop("No comparison results found; run compare_spatial_snrna_DE.R for each line first")
res <- lapply(names(found), function(l) plot_line(found[[l]], l))

## ---------------- Cross-line overview ----------------
all_summ <- bind_rows(lapply(res, `[[`, "summ")) %>%
  mutate(comparison = ifelse(snrna_group == "all", line, paste0(line, "\n", snrna_group)))
all_degs <- bind_rows(lapply(res, `[[`, "degs")) %>%
  mutate(comparison = ifelse(snrna_group == "all", line, paste0(line, "\n", snrna_group)))
lvl <- sort(unique(all_summ$comparison))

rho_stats <- all_summ %>% group_by(comparison) %>%
  summarise(n = n(), pos = sum(spearman_rho > 0), med = median(spearman_rho),
            p = binom.test(pos, n, 0.5, alternative = "greater")$p.value, .groups = "drop") %>%
  mutate(label = sprintf("%d/%d\np=%.2g", pos, n, p))

ov1 <- ggplot(all_summ, aes(factor(comparison, levels = lvl), spearman_rho)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.6) +
  geom_jitter(aes(colour = cell_class(spatial)), width = 0.15, size = 0.9, alpha = 0.8) +
  geom_text(data = rho_stats, aes(y = Inf, label = label), vjust = 1.2, size = 2.8) +
  labs(x = NULL, y = "Spearman rho per cell type", colour = "class",
       title = "Xenium vs snRNA-seq agreement across lines",
       subtitle = "one point per cell type; text = cell types with rho > 0, sign test p") +
  theme_bw(base_size = 10)
ggsave(file.path(overview_dir, "overview_rho_by_line.png"), ov1,
       width = max(6, 1.3 * length(lvl) + 2), height = 5, dpi = 300)

dir_stats <- all_degs %>% group_by(comparison) %>%
  summarise(n = n(), same = sum(same_direction),
            p = binom.test(same, n, 0.5, alternative = "greater")$p.value, .groups = "drop") %>%
  mutate(frac = same / n, label = sprintf("%d/%d\np=%.2g", same, n, p))

ov2 <- ggplot(dir_stats, aes(factor(comparison, levels = lvl), frac)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(aes(label = label), vjust = -0.3, size = 2.8) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "fraction of snRNA DEGs with same direction in Xenium",
       title = "Direction agreement of snRNA DEGs across lines",
       subtitle = "all cell types pooled; dashed = chance") +
  theme_bw(base_size = 10)
ggsave(file.path(overview_dir, "overview_direction_by_line.png"), ov2,
       width = max(6, 1.3 * length(lvl) + 2), height = 5, dpi = 300)

write_csv(full_join(rho_stats %>% select(comparison, n_celltypes = n, n_rho_positive = pos, median_rho = med, rho_sign_p = p),
                    dir_stats %>% select(comparison, n_snrna_degs = n, n_same_dir = same, frac_same_dir = frac, dir_sign_p = p),
                    by = "comparison"),
          file.path(overview_dir, "overview_stats.csv"))

pdf(file.path(overview_dir, "overview_figures.pdf"), width = 11, height = 8.5)
print(ov1); print(ov2)
dev.off()
message("Overview in ", overview_dir)