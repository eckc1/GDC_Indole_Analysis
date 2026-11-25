#!/usr/bin/env Rscript
# ==============================================================
# HALLMARK OVER-REPRESENTATION ANALYSIS (ORA)
# Indole + LPS vs LPS alone | DESeq2 significant genes only
#
# Test: one-sided hypergeometric enrichment (Fisher equivalent)
# Background: all eligible genes in the DESeq2 results that can
#             be mapped to gene symbols (not all human genes).
# Selection: DESeq2 padj < 0.05 (optional |log2FC| cutoff).
# Gene sets: MSigDB Human Hallmark (H) via msigdbr.
# Separate analyses: All significant, Up, and Down.
# Multiple testing: BH FDR within each analysis, across all
#                   Hallmark sets passing size/coverage filters.
# ==============================================================

# ---------- 0. User settings ----------
base_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole"
deg_file <- file.path(base_dir, "DESeq2_Paired_All20", "DEG_INLPS_vs_LPS_All.csv")
out_dir <- file.path(base_dir, "ORA")

DEG_FDR_CUTOFF <- 0.05
MIN_ABS_LOG2FC <- 0   # 0 = no additional fold-change threshold
ORA_FDR_CUTOFF <- 0.05
MIN_SET_SIZE <- 10    # size AFTER overlap with DESeq2 background
MAX_SET_SIZE <- 500
TOP_N_PLOT <- 20

# ---------- 1. Packages ----------
cran <- c("dplyr", "tidyr", "ggplot2", "msigdbr", "ggrepel")
for (pkg in cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

if (!file.exists(deg_file)) stop("DESeq2 input not found: ", deg_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 2. Load and standardize DESeq2 data ----------
deg <- read.csv(deg_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("log2FoldChange", "padj") %in% names(deg))) {
  stop("DESeq2 CSV needs log2FoldChange and padj columns")
}

symbol_field <- intersect(c("gene_name", "gene_symbol", "symbol", "SYMBOL"), names(deg))
if (length(symbol_field)) {
  deg$gene <- as.character(deg[[symbol_field[1]]])
  # Map rows missing symbols if gene_id is available.
  bad <- is.na(deg$gene) | !nzchar(trimws(deg$gene)) | toupper(trimws(deg$gene)) == "NA"
} else {
  if (!"gene_id" %in% names(deg)) stop("Provide gene symbol column or gene_id")
  deg$gene <- rep(NA_character_, nrow(deg))
  bad <- rep(TRUE, nrow(deg))
}
if (any(bad) && "gene_id" %in% names(deg)) {
  ids <- sub("\\.[0-9]+$", "", as.character(deg$gene_id[bad]))
  valid <- !is.na(ids) & nzchar(ids)
  mapped_symbols <- rep(NA_character_, length(ids))
  if (any(valid)) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(ids[valid]), column = "SYMBOL", keytype = "ENSEMBL",
      multiVals = "first"
    )
    mapped_symbols[valid] <- unname(mapped[ids[valid]])
  }
  deg$gene[bad] <- mapped_symbols
}

deg <- deg %>%
  mutate(
    gene = toupper(trimws(as.character(gene))),
    log2FC = suppressWarnings(as.numeric(log2FoldChange)),
    gene_FDR = suppressWarnings(as.numeric(padj)),
    gene_pvalue = if ("pvalue" %in% names(deg)) suppressWarnings(as.numeric(pvalue)) else NA_real_
  ) %>%
  filter(!is.na(gene), nzchar(gene), gene != "NA", is.finite(log2FC)) %>%
  # If several Ensembl IDs map to the same symbol, use the row
  # with the best finite adjusted p-value; preserve one background gene.
  arrange(is.na(gene_FDR), gene_FDR, is.na(gene_pvalue), gene_pvalue,
          desc(abs(log2FC)), gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FC, gene_FDR, gene_pvalue)

universe <- unique(deg$gene)
sig <- deg %>% filter(!is.na(gene_FDR), gene_FDR < DEG_FDR_CUTOFF,
                      abs(log2FC) >= MIN_ABS_LOG2FC)
gene_lists <- list(
  All_Significant = unique(sig$gene),
  Up_Indole_LPS = unique(sig$gene[sig$log2FC > 0]),
  Down_Indole_LPS = unique(sig$gene[sig$log2FC < 0])
)
if (!length(gene_lists$All_Significant)) stop("No significant DESeq2 genes at requested thresholds")
write.csv(deg, file.path(out_dir, "ORA_DESeq2_Background_Genes.csv"), row.names = FALSE)
write.csv(sig, file.path(out_dir, "ORA_DESeq2_Significant_Genes.csv"), row.names = FALSE)

# ---------- 3. Human Hallmark sets ----------
hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  transmute(Pathway = as.character(gs_name), gene = toupper(trimws(as.character(gene_symbol)))) %>%
  filter(!is.na(gene), nzchar(gene)) %>% distinct()
write.csv(hallmark_df, file.path(out_dir, "ORA_Hallmark_Gene_Set_Definitions.csv"), row.names = FALSE)
raw_sets <- split(hallmark_df$gene, hallmark_df$Pathway)
# Restrict both the pathway size and membership to the actual tested universe.
sets <- lapply(raw_sets, function(x) intersect(unique(x), universe))
set_sizes <- lengths(sets)
sets <- sets[set_sizes >= MIN_SET_SIZE & set_sizes <= MAX_SET_SIZE]
if (!length(sets)) stop("No Hallmark sets meet the size filters within the DESeq2 universe")

# ---------- 4. One-sided hypergeometric ORA ----------
# N = tested background, K = genes in set, n = selected significant genes,
# k = selected significant genes in set.
# P(X >= k) = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
run_ora <- function(selected_genes, analysis_name) {
  selected_genes <- intersect(unique(selected_genes), universe)
  N <- length(universe)
  n <- length(selected_genes)
  rows <- lapply(names(sets), function(pathway) {
    members <- sets[[pathway]]
    K <- length(members)
    overlap <- intersect(selected_genes, members)
    k <- length(overlap)
    p <- if (n == 0 || k == 0) 1 else
      phyper(k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    expected <- n * K / N
    data.frame(
      Analysis = analysis_name, Pathway = pathway,
      Background_Size = N, Selected_Genes = n,
      Pathway_Size_In_Background = K, Overlap_Count = k,
      Expected_Overlap = expected,
      Fold_Enrichment = if (expected > 0) k / expected else NA_real_,
      P_Value = p,
      Overlap_Genes = paste(sort(overlap), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  res <- dplyr::bind_rows(rows)
  res$FDR <- p.adjust(res$P_Value, method = "BH")
  res$Significant <- res$FDR < ORA_FDR_CUTOFF
  res$Pathway_Label <- gsub("_", " ", sub("^HALLMARK_", "", res$Pathway))
  res %>% arrange(FDR, P_Value, desc(Fold_Enrichment))
}

results <- dplyr::bind_rows(lapply(names(gene_lists), function(group) {
  run_ora(gene_lists[[group]], group)
}))
write.csv(results, file.path(out_dir, "ORA_Hallmark_All_Results.csv"), row.names = FALSE)
write.csv(results %>% filter(Significant),
          file.path(out_dir, "ORA_Hallmark_Significant_Results.csv"), row.names = FALSE)

# Flatten each pathway gene overlap to an independently filterable gene-level file.
gene_detail <- results %>%
  filter(Overlap_Count > 0) %>%
  dplyr::select(Analysis, Pathway, FDR, Fold_Enrichment, Overlap_Genes) %>%
  separate_rows(Overlap_Genes, sep = ",\\s*") %>%
  rename(gene = Overlap_Genes) %>%
  left_join(deg, by = "gene") %>%
  arrange(Analysis, FDR, Pathway, gene_FDR)
write.csv(gene_detail, file.path(out_dir, "ORA_Hallmark_Overlapping_Genes_DESeq2.csv"),
          row.names = FALSE)

# ---------- 5. Figures (per analysis) ----------
plot_dir <- file.path(out_dir, "Figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

for (group in names(gene_lists)) {
  this <- results %>% filter(Analysis == group, Overlap_Count > 0)
  if (!nrow(this)) {
    message("No overlapping genes to plot for ", group)
    next
  }
  significant <- this %>% filter(Significant)
  to_plot <- if (nrow(significant)) significant else this
  plot_note <- if (nrow(significant)) {
    sprintf("Hallmark ORA; BH FDR < %.2f", ORA_FDR_CUTOFF)
  } else {
    "Exploratory: no Hallmark pathways passed the ORA FDR cutoff"
  }
  to_plot <- to_plot %>% arrange(FDR, P_Value) %>% slice_head(n = TOP_N_PLOT)
  to_plot <- to_plot %>% mutate(
    Plot_Label = factor(Pathway_Label, levels = rev(Pathway_Label)),
    MinusLog10FDR = -log10(pmax(FDR, .Machine$double.xmin))
  )
  h <- max(5.5, 0.38 * nrow(to_plot) + 2.0)

  # Figure 1: fold-enrichment bars (unlike the GSEA gene-segment bars,
  # these bar lengths are statistically interpretable enrichment ratios).
  p_bar <- ggplot(to_plot, aes(x = Fold_Enrichment, y = Plot_Label)) +
    geom_col(fill = "#4477AA", width = 0.75) +
    geom_text(aes(label = sprintf("%d genes | q %.2g", Overlap_Count, FDR)), hjust = -0.08, size = 3.1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(title = paste("Hallmark ORA -", gsub("_", " ", group)),
         subtitle = plot_note, x = "Fold enrichment (observed / expected)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(plot_dir, paste0("01_ORA_Fold_Enrichment_", group, ".png")),
         p_bar, width = 12, height = h, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0("01_ORA_Fold_Enrichment_", group, ".pdf")),
         p_bar, width = 12, height = h)

  # Figure 2: significance + gene-overlap dotplot.
  p_dot <- ggplot(to_plot, aes(x = Fold_Enrichment, y = Plot_Label)) +
    geom_point(aes(size = Overlap_Count, color = MinusLog10FDR), alpha = 0.9) +
    scale_color_gradient(low = "#4477AA", high = "#BB3322") +
    labs(title = paste("Hallmark ORA -", gsub("_", " ", group)),
         subtitle = plot_note, x = "Fold enrichment (observed / expected)",
         y = NULL, size = "Overlapping DEGs", color = "-log10(ORA FDR)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(plot_dir, paste0("02_ORA_Dotplot_", group, ".png")),
         p_dot, width = 12, height = h, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0("02_ORA_Dotplot_", group, ".pdf")),
         p_dot, width = 12, height = h)
}


# ---------- 6. Combined ORA figures and ONE GLOBAL VOLCANO per cutoff ----------
# DESeq2 gene threshold stays at 0.05 for BOTH analyses.
# Only ORA pathway FDR changes: 0.05 versus 0.10.
# No individual pathway volcano plots are generated.
ORA_PLOT_CUTOFFS <- c(0.05, 0.10)
GENES_PER_BAR <- 10
GENES_PER_VOLCANO <- 5

volcano_bg <- deg %>%
  mutate(volcano_p = if_else(is.finite(gene_pvalue) & gene_pvalue > 0,
                             gene_pvalue, gene_FDR),
         volcano_p = pmax(volcano_p, .Machine$double.xmin),
         minus_log10_p = -log10(volcano_p)) %>%
  filter(is.finite(minus_log10_p))

make_ora_volcano <- function(highlight_genes, plot_title, plot_subtitle) {
  hits <- volcano_bg %>%
    semi_join(highlight_genes %>% distinct(gene), by = "gene") %>%
    mutate(Highlight_Direction = if_else(log2FC > 0,
                                       "Higher in Indole + LPS", "Higher in LPS"))
  ggplot(volcano_bg, aes(x = log2FC, y = minus_log10_p)) +
    geom_point(color = "grey75", size = 0.7, alpha = 0.45) +
    geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
    geom_point(data = hits, aes(color = Highlight_Direction),
               size = 2.2, alpha = 0.95) +
    ggrepel::geom_text_repel(data = hits, aes(label = gene),
                            size = 3, max.overlaps = Inf,
                            max.time = 5, max.iter = 20000,
                            box.padding = 0.32, point.padding = 0.18,
                            segment.color = NA, seed = 123) +
    scale_color_manual(values = c("Higher in Indole + LPS" = "#D73027",
                                  "Higher in LPS" = "#4575B4")) +
    labs(title = plot_title, subtitle = plot_subtitle,
         x = "DESeq2 log2 fold change (Indole + LPS / LPS)",
         y = "-log10(DESeq2 raw p-value; FDR fallback)",
         color = "Highlighted gene") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
}

dir_colors <- c("Higher in Indole + LPS" = "#D73027",
                "Higher in LPS" = "#4575B4")

for (plot_cutoff in ORA_PLOT_CUTOFFS) {
  tag <- sprintf("FDR_%.2f", plot_cutoff)
  cutoff_dir <- file.path(plot_dir, tag)
  dir.create(cutoff_dir, recursive = TRUE, showWarnings = FALSE)
  message("Generating combined ORA charts and global volcano for ", tag)

  combined <- results %>%
    filter(Analysis %in% c("Up_Indole_LPS", "Down_Indole_LPS"),
           Overlap_Count > 0, is.finite(Fold_Enrichment),
           !is.na(FDR), FDR < plot_cutoff) %>%
    group_by(Analysis) %>%
    arrange(FDR, P_Value, .by_group = TRUE) %>%
    slice_head(n = TOP_N_PLOT) %>% ungroup() %>%
    mutate(Direction = if_else(Analysis == "Up_Indole_LPS",
                               "Higher in Indole + LPS", "Higher in LPS"),
           Signed_Enrichment = if_else(Analysis == "Up_Indole_LPS",
                                       Fold_Enrichment, -Fold_Enrichment),
           MinusLog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
           Pathway_Direction = paste(Pathway_Label, Direction, sep = " | ")) %>%
    arrange(Signed_Enrichment) %>%
    mutate(Pathway_Direction = factor(Pathway_Direction,
                                      levels = unique(Pathway_Direction)))

  write.csv(combined, file.path(cutoff_dir, "ORA_Combined_Up_Down_Pathways.csv"),
            row.names = FALSE)

  if (nrow(combined)) {
    h <- max(7, 0.33 * nrow(combined) + 2.2)
    p_bar <- ggplot(combined,
                    aes(x = Signed_Enrichment, y = Pathway_Direction)) +
      geom_vline(xintercept = 0, color = "grey45", linewidth = 0.4) +
      geom_col(aes(fill = Direction), width = 0.72) +
     # geom_text(aes(label = sprintf("%d genes | q %.2g", Overlap_Count, FDR),
                  #  hjust = if_else(Signed_Enrichment > 0, -0.08, 1.08)),
                #size = 2.8) +
      scale_fill_manual(values = dir_colors) +
      scale_x_continuous(labels = function(x) format(abs(x), trim = TRUE),
                         expand = expansion(mult = c(0.24, 0.24))) +
      labs(title = "Hallmark ORA: Upregulated vs Downregulated DEGs",
           subtitle = sprintf("DESeq2 FDR < %.2f | ORA FDR < %.2f",
                              DEG_FDR_CUTOFF, plot_cutoff),
           x = "Fold enrichment (left: LPS higher | right: Indole + LPS higher)",
           y = NULL, fill = "Gene direction") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(size = 9), legend.position = "bottom")
    ggsave(file.path(cutoff_dir, "03_ORA_Combined_Up_Down_Barplot.png"),
           p_bar, width = 17, height = h, dpi = 300, bg = "white", limitsize = FALSE)
    ggsave(file.path(cutoff_dir, "03_ORA_Combined_Up_Down_Barplot.pdf"),
           p_bar, width = 17, height = h, limitsize = FALSE)

    p_dot <- ggplot(combined,
                    aes(x = Signed_Enrichment, y = Pathway_Direction)) +
      geom_vline(xintercept = 0, color = "grey45", linewidth = 0.4) +
      geom_point(aes(color = Direction, size = Overlap_Count,
                     alpha = MinusLog10FDR)) +
      scale_color_manual(values = dir_colors) +
      scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
      scale_x_continuous(labels = function(x) format(abs(x), trim = TRUE),
                         expand = expansion(mult = c(0.12, 0.12))) +
      labs(title = "Hallmark ORA: Upregulated vs Downregulated DEGs",
           subtitle = sprintf("ORA FDR < %.2f", plot_cutoff),
           x = "Fold enrichment (left: LPS higher | right: Indole + LPS higher)",
           y = NULL, color = "Gene direction", size = "Overlapping DEGs") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(size = 9), legend.position = "bottom")
    ggsave(file.path(cutoff_dir, "04_ORA_Combined_Up_Down_Dotplot.png"),
           p_dot, width = 17, height = h, dpi = 300, bg = "white", limitsize = FALSE)
    ggsave(file.path(cutoff_dir, "04_ORA_Combined_Up_Down_Dotplot.pdf"),
           p_dot, width = 17, height = h, limitsize = FALSE)
  }

  gene_ranked <- gene_detail %>%
    filter(Analysis %in% c("Up_Indole_LPS", "Down_Indole_LPS")) %>%
    semi_join(combined %>% dplyr::select(Analysis, Pathway), by = c("Analysis", "Pathway")) %>%
    filter(!is.na(gene_FDR), gene_FDR < DEG_FDR_CUTOFF, is.finite(log2FC)) %>%
    group_by(Analysis, Pathway) %>%
    arrange(gene_FDR, is.na(gene_pvalue), gene_pvalue,
            desc(abs(log2FC)), gene, .by_group = TRUE) %>%
    mutate(Gene_Rank = row_number()) %>% ungroup() %>%
    left_join(combined %>%
                dplyr::select(Analysis, Pathway, Direction, ORA_FDR = FDR,
                       ORA_Fold_Enrichment = Fold_Enrichment),
              by = c("Analysis", "Pathway"))

  top10 <- gene_ranked %>% filter(Gene_Rank <= GENES_PER_BAR)
  top5 <- gene_ranked %>% filter(Gene_Rank <= GENES_PER_VOLCANO)
  write.csv(gene_ranked, file.path(cutoff_dir, "ORA_Selected_Pathways_All_Significant_Genes.csv"), row.names = FALSE)
  write.csv(top10, file.path(cutoff_dir, "ORA_Top10_Genes_Per_Pathway.csv"), row.names = FALSE)
  write.csv(top5, file.path(cutoff_dir, "ORA_Top5_Genes_Per_Pathway_Volcano.csv"), row.names = FALSE)

  if (nrow(top10)) {
    bar_labels <- combined %>%
      mutate(Display_Label = sprintf("%s | %s | ORA q %.2g",
                                     Pathway_Label,
                                     if_else(Analysis == "Up_Indole_LPS", "UP", "DOWN"), FDR)) %>%
      dplyr::select(Analysis, Pathway, Display_Label, Signed_Enrichment) %>%
      arrange(Signed_Enrichment)
    bars <- top10 %>%
      left_join(bar_labels %>% dplyr::select(Analysis, Pathway, Display_Label),
                by = c("Analysis", "Pathway")) %>%
      mutate(Display_Label = factor(Display_Label, levels = bar_labels$Display_Label))
    p_genes <- ggplot(bars, aes(x = Gene_Rank, y = Display_Label)) +
      geom_tile(aes(fill = Direction), width = 0.95, height = 0.83,
                color = "white", linewidth = 0.35) +
      geom_text(aes(label = gene), color = "white", fontface = "bold", size = 2.55) +
      scale_fill_manual(values = dir_colors) +
      scale_x_continuous(breaks = seq_len(GENES_PER_BAR),
                         limits = c(0.5, GENES_PER_BAR + 0.5),
                         expand = expansion(mult = 0)) +
      labs(title = "Hallmark ORA: Top 10 significant DESeq2 genes per pathway",
           subtitle = sprintf("DESeq2 FDR < %.2f | ORA FDR < %.2f | Ranked by DESeq2 FDR",
                              DEG_FDR_CUTOFF, plot_cutoff),
           x = "Gene rank (1 = lowest DESeq2 adjusted p-value)",
           y = NULL, fill = "Gene direction") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9),
            legend.position = "bottom", plot.title = element_text(face = "bold"))
    bh <- max(7, 0.45 * nrow(bar_labels) + 2)
    ggsave(file.path(cutoff_dir, "05_ORA_Top10_Gene_Bars_Combined.png"),
           p_genes, width = 19, height = bh, dpi = 300, limitsize = FALSE, bg = "white")
    ggsave(file.path(cutoff_dir, "05_ORA_Top10_Gene_Bars_Combined.pdf"),
           p_genes, width = 19, height = bh, limitsize = FALSE)
  }

  # ONE global volcano per cutoff. Each unique selected gene appears once,
  # even if it belongs to multiple ORA-enriched pathways.
  global_genes <- top5 %>% distinct(gene)
  write.csv(global_genes, file.path(cutoff_dir, "ORA_Global_Top5_Unique_Genes.csv"), row.names = FALSE)
  global_plot <- make_ora_volcano(
    global_genes,
    sprintf("Indole + LPS vs LPS: ORA top-5 genes (FDR < %.2f)", plot_cutoff),
    sprintf("%d unique highlighted genes across %d pathway-direction pairs",
            nrow(global_genes), nrow(combined)))
  ggsave(file.path(cutoff_dir, "06_ORA_Global_Volcano_Top5.png"),
         global_plot, width = 15, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(cutoff_dir, "06_ORA_Global_Volcano_Top5.pdf"),
         global_plot, width = 15, height = 11)
  message(tag, ": ", nrow(combined), " pathway-direction pairs, ",
          nrow(global_genes), " unique highlighted genes. Saved to ", cutoff_dir)
}

# ---------- 7. Session summary ----------
capture.output({
  cat("Human Hallmark ORA: Indole + LPS vs LPS\n")
  cat("DESeq2 gene cutoff:", DEG_FDR_CUTOFF, "\n")
  cat("ORA plotting cutoffs:", paste(ORA_PLOT_CUTOFFS, collapse = ", "), "\n")
  cat("Background genes:", length(universe), "\n")
  cat("MSigDB Hallmark sets:", length(sets), "\n")
  print(sessionInfo())
}, file = file.path(out_dir, "ORA_Analysis_Summary_Dual_FDR.txt"))
message("ORA COMPLETE. Dual-cutoff charts: ", plot_dir)
