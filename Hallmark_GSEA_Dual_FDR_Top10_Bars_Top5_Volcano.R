#!/usr/bin/env Rscript
# ============================================================
# HALLMARK GSEA + DESeq2: DUAL PATHWAY FDR (0.05/0.10): TOP 10 BAR GENES + TOP 5 VOLCANO GENES
#
# Inputs (from your existing scripts):
#   FDR_0.1_Hallmark_GSEA/Hallmark_GSEA_All_Results.csv
#   DESeq2_Paired_All20/DEG_INLPS_vs_LPS_All.csv
#
# Outputs:
#   Top10_Genes_Per_Pathway_All.csv
#   Top10_Genes_Per_Pathway_Plotted.csv
#   All_Leading_Edge_Genes_Merged_DESeq2.csv
#   01_Pathway_Top10_Gene_Bars.png/.pdf
#   02_Global_Volcano_Selected_Top10.png/.pdf
#   03_Pathway_Volcanoes/*.png/.pdf
#
# Each colored, 10-segment horizontal bar represents the selected
# genes (NOT NES magnitude). Pathway NES/FDR are displayed beside it.
# Volcano plots use all tested DESeq2 genes as background, highlighting
# only the selected top genes for each pathway.
# ============================================================

# ---- 0. Settings ----
base_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole"
deg_file <- file.path(base_dir, "DESeq2_Paired_All20", "DEG_INLPS_vs_LPS_All.csv")
# Use one UNFILTERED GSEA results CSV containing all tested pathways.
# The 0.05 and 0.10 cutoffs are applied separately below.
gsea_file <- file.path(base_dir, "FDR_0.1_Hallmark_GSEA", "Hallmark_GSEA_All_Results.csv")
output_root <- file.path(base_dir, "FDR_0.1_Hallmark_GSEA", "Dual_FDR_Gene_Annotations")

pathway_fdr_cutoffs <- c(0.05, 0.10)
bar_n_genes <- 10
volcano_n_genes <- 5
max_pathways_per_direction <- 15
# TRUE: use leading-edge genes; FALSE: all Hallmark members.
use_leading_edge <- TRUE

# ---- 1. Dependencies ----
cran <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "readr", "stringr")
for (pkg in cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
if (!use_leading_edge && !requireNamespace("msigdbr", quietly = TRUE))
  install.packages("msigdbr", repos = "https://cloud.r-project.org")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)

stopifnot(file.exists(deg_file), file.exists(gsea_file))

# ---- 2. Read & standardize DESeq2 ----
deg <- read.csv(deg_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("log2FoldChange", "padj") %in% names(deg)))
  stop("DESeq2 CSV must have log2FoldChange and padj columns")

symbol_candidates <- c("gene_name", "gene_symbol", "symbol", "SYMBOL")
symbol_field <- intersect(symbol_candidates, names(deg))
if (length(symbol_field)) {
  deg$gene <- as.character(deg[[symbol_field[[1]]]])
} else {
  if (!"gene_id" %in% names(deg)) stop("Need a gene symbol column or gene_id")
  ids <- sub("\\.[0-9]+$", "", as.character(deg$gene_id))
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids),
    keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  deg$gene <- unname(mapped[ids])
}

# Keep one unique entry per gene symbol: best adjusted p value first.
# This avoids one-to-many symbol mapping duplicating pathway hits.
deg <- deg %>%
  mutate(gene = toupper(trimws(gene)),
         log2FC = suppressWarnings(as.numeric(log2FoldChange)),
         gene_FDR = suppressWarnings(as.numeric(padj)),
         gene_pvalue = if ("pvalue" %in% names(.))
           suppressWarnings(as.numeric(pvalue)) else NA_real_,
         gene_stat = if ("stat" %in% names(.))
           suppressWarnings(as.numeric(stat)) else NA_real_) %>%
  filter(!is.na(gene), nzchar(gene), gene != "NA", is.finite(log2FC)) %>%
  arrange(is.na(gene_FDR), gene_FDR, is.na(gene_pvalue), gene_pvalue,
          desc(abs(log2FC))) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FC, gene_FDR, gene_pvalue, gene_stat)

# ---- 3. Read GSEA pathways and get their genes ----
gsea <- read.csv(gsea_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("Pathway", "NES", "FDR") %in% names(gsea)))
  stop("Expected Pathway, NES, FDR in GSEA CSV")

gsea <- gsea %>%
  mutate(NES = as.numeric(NES), FDR = as.numeric(FDR),
         Direction = if_else(NES > 0, "Indole + LPS", "LPS"),
         Pathway_Label = str_replace_all(str_remove(Pathway, "^HALLMARK_"), "_", " ")) %>%
  filter(!is.na(Pathway), is.finite(NES)) %>%
  distinct(Pathway, .keep_all = TRUE)

if (use_leading_edge) {
  if (!"Leading_Edge_Genes" %in% names(gsea))
    stop("GSEA CSV lacks Leading_Edge_Genes. Set use_leading_edge = FALSE.")
  members <- gsea %>%
    dplyr::select(Pathway, Leading_Edge_Genes) %>%
    mutate(Leading_Edge_Genes = if_else(is.na(Leading_Edge_Genes),
                                        "", Leading_Edge_Genes)) %>%
    separate_rows(Leading_Edge_Genes, sep = ",\\s*") %>%
    transmute(Pathway, gene = toupper(trimws(Leading_Edge_Genes))) %>%
    filter(nzchar(gene), gene != "NA") %>% distinct()
} else {
  hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  members <- hallmark %>%
    transmute(Pathway = gs_name, gene = toupper(gene_symbol)) %>% distinct()
}

# This table retains ALL matched leading-edge genes, not just the top 10.
merged <- members %>%
  inner_join(gsea %>% dplyr::select(Pathway, Pathway_Label, NES, FDR, Direction),
             by = "Pathway") %>%
  inner_join(deg, by = "gene") %>%
  mutate(gene_significant = !is.na(gene_FDR) & gene_FDR < 0.05)

if (!nrow(merged)) stop("No gene matches. Check gene IDs and GSEA source data.")


# ---- 4. Sort leading-edge genes for each pathway ----
# Rank by gene-level adjusted p-value, breaking ties by raw p-value,
# then absolute log2 fold change. Missing adjusted p-values are last.
ranked_genes <- merged %>%
  group_by(Pathway) %>%
  arrange(is.na(gene_FDR), gene_FDR,
          is.na(gene_pvalue), gene_pvalue,
          desc(abs(log2FC)), .by_group = TRUE) %>%
  mutate(gene_rank_in_pathway = row_number()) %>%
  ungroup()

# ---- 5. Volcano background: all DESeq2 genes ----
volcano <- deg %>%
  mutate(volcano_p = if_else(is.finite(gene_pvalue) & gene_pvalue > 0,
                              gene_pvalue, gene_FDR),
         volcano_p = pmax(volcano_p, .Machine$double.xmin),
         minus_log10_p = -log10(volcano_p)) %>%
  filter(is.finite(minus_log10_p))

make_volcano <- function(highlights, plot_title, plot_subtitle) {
  hits <- volcano %>% semi_join(highlights %>% distinct(gene), by = "gene") %>%
    mutate(highlight_direction = if_else(log2FC >= 0,
                                         "Higher in Indole + LPS", "Higher in LPS"))
  ggplot(volcano, aes(log2FC, minus_log10_p)) +
    geom_point(color = "grey75", size = 0.7, alpha = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(data = hits, aes(color = highlight_direction),
               size = 2.1, alpha = 0.95) +
    ggrepel::geom_text_repel(data = hits, aes(label = gene),
                             size = 3.1, max.overlaps = Inf,
                             max.time = 5, max.iter = 20000,
                             box.padding = 0.3, point.padding = 0.18,
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

# ---- 6. Run both pathway-level FDR cutoffs ----
for (pathway_fdr_cutoff in pathway_fdr_cutoffs) {
  cutoff_label <- sprintf("FDR_%0.2f", pathway_fdr_cutoff)
  out_dir <- file.path(output_root, cutoff_label)
  volcano_dir <- file.path(out_dir, "03_Pathway_Volcanoes")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

  # Save merged genes and per-pathway rankings for auditing.
  write.csv(merged, file.path(out_dir, "All_Leading_Edge_Genes_Merged_DESeq2.csv"),
            row.names = FALSE)
  top10 <- ranked_genes %>% filter(gene_rank_in_pathway <= bar_n_genes)
  top5 <- ranked_genes %>% filter(gene_rank_in_pathway <= volcano_n_genes)
  write.csv(top10, file.path(out_dir, "Top10_Genes_Per_Pathway_All.csv"), row.names = FALSE)
  write.csv(top5, file.path(out_dir, "Top5_Genes_Per_Pathway_All.csv"), row.names = FALSE)

  eligible <- gsea %>% filter(!is.na(FDR), FDR < pathway_fdr_cutoff)
  if (!nrow(eligible)) {
    warning("No significant pathways at ", cutoff_label, "; output tables saved, no figures.")
    next
  }
  selected <- eligible %>%
    group_by(Direction) %>%
    arrange(FDR, desc(abs(NES)), .by_group = TRUE) %>%
    slice_head(n = max_pathways_per_direction) %>% ungroup()

  plot_top10 <- top10 %>% semi_join(selected, by = "Pathway")
  plot_top5 <- top5 %>% semi_join(selected, by = "Pathway")
  if (!nrow(plot_top10) || !nrow(plot_top5)) {
    warning("No matching DESeq2 genes for pathways at ", cutoff_label)
    next
  }
  write.csv(plot_top10, file.path(out_dir, "Top10_Genes_Per_Pathway_Plotted.csv"),
            row.names = FALSE)
  write.csv(plot_top5, file.path(out_dir, "Top5_Genes_Per_Pathway_Plotted.csv"),
            row.names = FALSE)
  write.csv(selected, file.path(out_dir, "Selected_Pathways.csv"), row.names = FALSE)

  # ---- 7. Segmented pathway bar chart: top TEN genes ----
  pathway_labels <- plot_top10 %>%
    distinct(Pathway, Pathway_Label, NES, FDR) %>%
    arrange(NES) %>%
    mutate(label = sprintf("%s   |   NES %+.2f   |   q %.2g",
                           Pathway_Label, NES, FDR))
  plot_bar <- plot_top10 %>%
    left_join(pathway_labels %>% dplyr::select(Pathway, label), by = "Pathway") %>%
    mutate(label = factor(label, levels = pathway_labels$label),
           direction = factor(Direction, levels = c("LPS", "Indole + LPS")))

  p_bar <- ggplot(plot_bar, aes(x = gene_rank_in_pathway, y = label)) +
    geom_tile(aes(fill = direction), width = 0.96, height = 0.84,
              color = "white", linewidth = 0.5) +
    geom_text(aes(label = gene), size = 2.7, color = "white", fontface = "bold") +
    scale_x_continuous(breaks = seq_len(bar_n_genes),
                       limits = c(0.5, bar_n_genes + 0.5),
                       expand = expansion(mult = 0)) +
    scale_fill_manual(values = c("Indole + LPS" = "#D73027", "LPS" = "#4575B4")) +
    labs(title = "Top 10 DESeq2 genes within each enriched Hallmark pathway",
         subtitle = sprintf("Hallmark pathways FDR < %.2f | Leading-edge genes ranked by DESeq2 adjusted p-value", pathway_fdr_cutoff),
         x = "Gene rank within pathway (1 = smallest gene-level adjusted p-value)",
         y = NULL, fill = "GSEA enrichment toward") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9),
          legend.position = "bottom", plot.title = element_text(face = "bold"))

  bar_height <- max(7, 0.43 * n_distinct(plot_bar$Pathway) + 2)
  ggsave(file.path(out_dir, "01_Pathway_Top10_Gene_Bars.png"), p_bar,
         width = 19, height = bar_height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(file.path(out_dir, "01_Pathway_Top10_Gene_Bars.pdf"), p_bar,
         width = 19, height = bar_height, limitsize = FALSE)

  # ---- 8. Global volcano: only top FIVE genes per selected pathway ----
  global_hits <- plot_top5 %>% distinct(gene)
  p_global <- make_volcano(
    global_hits,
    "Indole + LPS vs LPS: GSEA top-5 leading-edge genes",
    sprintf("GSEA FDR < %.2f | %d unique genes across %d pathways",
            pathway_fdr_cutoff, nrow(global_hits), n_distinct(plot_top5$Pathway)))
  ggsave(file.path(out_dir, "02_Global_Volcano_Selected_Top5.png"), p_global,
         width = 15, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "02_Global_Volcano_Selected_Top5.pdf"), p_global,
         width = 15, height = 11)

  # ---- 9. One volcano per pathway: top FIVE genes ----
  for (pathway_name in unique(plot_top5$Pathway)) {
    this_pathway <- plot_top5 %>% filter(Pathway == pathway_name)
    p <- make_volcano(
      this_pathway,
      paste0("Volcano: ", unique(this_pathway$Pathway_Label)),
      sprintf("GSEA cutoff %.2f | NES %+.3f | pathway FDR %.3g | %d highlighted genes",
              pathway_fdr_cutoff, unique(this_pathway$NES),
              unique(this_pathway$FDR), nrow(this_pathway)))
    safe <- gsub("[^A-Za-z0-9_]+", "_", pathway_name)
    ggsave(file.path(volcano_dir, paste0(safe, "_Top5_Volcano.png")), p,
           width = 11, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(volcano_dir, paste0(safe, "_Top5_Volcano.pdf")), p,
           width = 11, height = 8)
  }
  message("COMPLETE ", cutoff_label, ": ", n_distinct(plot_top10$Pathway),
          " plotted pathways, ", nrow(global_hits), " unique top-5 volcano genes. Output: ", out_dir)
}
