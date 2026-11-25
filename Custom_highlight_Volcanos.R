
#!/usr/bin/env Rscript

# ============================================================
# HUMAN GDC INDOLE + LPS vs LPS
# TARGETED VOLCANO PLOT
#
# Highlights:
#   1. AhR signaling and response genes
#   2. Cytokines and chemokines
#   3. Immune activation markers
#
# Rules:
#   - Only highlight genes with FDR < 0.05
#   - Maximum 5 UP genes per category
#   - Maximum 5 DOWN genes per category
#   - Rank genes by smallest FDR
#   - All other genes appear gray
#   - Only highlighted genes receive labels
#
# Maximum highlighted genes = 30
#
# Output:
#   Targeted_Volcano_AhR_Cytokines_Activation.png
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

packages <- c(
  "dplyr",
  "ggplot2",
  "ggrepel"
)

missing <- packages[
  !vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

library(dplyr)
library(ggplot2)
library(ggrepel)


# ============================================================
# 1. INPUT AND OUTPUT
# ============================================================

base_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole"

input_file <- file.path(
  base_dir,
  "DESeq2_Paired_All20",
  "DEG_INLPS_vs_LPS_All.csv"
)

output_dir <- file.path(
  base_dir,
  "Targeted_Volcano_Plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_dir,
  "Targeted_Volcano_AhR_Cytokines_Activation.png"
)


# ============================================================
# 2. SETTINGS
# ============================================================

fdr_cutoff <- 0.05

# Maximum selected genes per category per direction
max_genes_per_direction <- 10

# Minimum fold change for selecting genes
# 0 means any nonzero fold change is allowed.
min_abs_log2fc <- 0

plot_width <- 15
plot_height <- 12

plot_dpi <- 300


# ============================================================
# 3. DEFINE TARGET GENE SETS
# ============================================================

# ------------------------------------------------------------
# A. AhR SIGNALING AND RESPONSE GENES
# ------------------------------------------------------------

ahr_genes <- c(
  
  # Core AhR signaling
  "AHR",
  "ARNT",
  "AHRR",
  
  # Canonical xenobiotic response
  "CYP1A1",
  "CYP1A2",
  "CYP1B1",
  "TIPARP",
  "NQO1",
  "ALDH3A1",
  
  # Receptor regulation
  "AIP",
  "HSP90AA1",
  "HSP90AB1",
  "PTGES3",
  
  # Additional AhR-associated response genes
  "UGT1A1",
  "GSTA1",
  "GSTM1",
  "GSTM2",
  "GSTP1",
  "IDO1",
  "TDO2"
)


# ------------------------------------------------------------
# B. CYTOKINES AND CHEMOKINES
# ------------------------------------------------------------

cytokine_genes <- c(
  
  # Pro-inflammatory cytokines
  "IL1A",
  "IL1B",
  "IL6",
  "TNF",
  "IL18",
  "IL12A",
  "IL12B",
  "IL23A",
  
  # Regulatory cytokines
  "IL10",
  "TGFB1",
  "IL1RN",
  
  # Interferons
  "IFNB1",
  "IFNG",
  "IFNA1",
  
  # CXCL chemokines
  "CXCL1",
  "CXCL2",
  "CXCL3",
  "CXCL5",
  "CXCL8",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CXCL12",
  
  # CCL chemokines
  "CCL2",
  "CCL3",
  "CCL4",
  "CCL5",
  "CCL7",
  "CCL8",
  "CCL20",
  
  # Colony stimulating factors
  "CSF1",
  "CSF2",
  "CSF3",
  
  # Other cytokines
  "LIF",
  "OSM",
  "IL15",
  "IL16",
  "IL27"
)


# ------------------------------------------------------------
# C. IMMUNE ACTIVATION MARKERS
# ------------------------------------------------------------

activation_genes <- c(
  
  # Costimulatory markers
  "CD80",
  "CD86",
  "CD40",
  "CD83",
  "CD69",
  "CD70",
  
  # Antigen presentation
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQB1",
  "CD74",
  
  # Activation / adhesion
  "ICAM1",
  "VCAM1",
  "CD274",
  "FAS",
  "CCR7",
  
  # Inflammatory signaling
  "NFKBIA",
  "NFKB1",
  "RELA",
  "RELB",
  "TNFAIP3",
  "STAT1",
  "STAT3",
  "SOCS1",
  "SOCS3",
  "IRF1",
  "IRF7",
  
  # Myeloid activation
  "CD14",
  "FCGR1A",
  "FCGR3A",
  "TLR2",
  "TLR4",
  "CLEC4E",
  "CD163",
  "CD68"
)


# ============================================================
# 4. BUILD GENE REFERENCE TABLE
# ============================================================

gene_reference <- bind_rows(
  
  data.frame(
    gene = ahr_genes,
    Category = "AhR Signaling"
  ),
  
  data.frame(
    gene = cytokine_genes,
    Category = "Cytokines"
  ),
  
  data.frame(
    gene = activation_genes,
    Category = "Activation Markers"
  )
  
) %>%
  mutate(
    gene = toupper(trimws(gene))
  ) %>%
  distinct(gene, .keep_all = TRUE)

cat(
  "\nTarget genes defined:",
  nrow(gene_reference),
  "\n"
)


# ============================================================
# 5. READ DESEQ2 RESULTS
# ============================================================

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

deg <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "log2FoldChange",
  "padj"
)

if (!all(required_columns %in% names(deg))) {
  stop("Missing log2FoldChange or padj column.")
}


# ============================================================
# 6. GENE SYMBOL ANNOTATION
# ============================================================

# Use existing gene symbols when available.
# Otherwise convert Ensembl IDs to symbols.

symbol_columns <- intersect(
  c(
    "gene_name",
    "gene_symbol",
    "symbol",
    "SYMBOL"
  ),
  names(deg)
)

if (length(symbol_columns) > 0) {
  
  deg$gene <- as.character(
    deg[[symbol_columns[1]]]
  )
  
} else {
  
  if (!"gene_id" %in% names(deg)) {
    stop("No gene symbol or gene_id column found.")
  }
  
  if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
  )) {
    install.packages("BiocManager")
  }
  
  bioc_packages <- c(
    "AnnotationDbi",
    "org.Hs.eg.db"
  )
  
  for (pkg in bioc_packages) {
    
    if (!requireNamespace(
      pkg,
      quietly = TRUE
    )) {
      
      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      )
      
    }
  }
  
  ensembl_ids <- sub(
    "\\.[0-9]+$",
    "",
    as.character(deg$gene_id)
  )
  
  mapping <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ensembl_ids),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  deg$gene <- unname(
    mapping[ensembl_ids]
  )
}

deg$gene <- toupper(
  trimws(as.character(deg$gene))
)


# ============================================================
# 7. CLEAN AND PREPARE DE RESULTS
# ============================================================

deg <- deg %>%
  mutate(
    log2FC = as.numeric(log2FoldChange),
    FDR = as.numeric(padj)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "NA",
    is.finite(log2FC)
  ) %>%
  arrange(
    is.na(FDR),
    FDR
  ) %>%
  distinct(
    gene,
    .keep_all = TRUE
  ) %>%
  mutate(
    
    # Avoid infinite -log10 values for FDR = 0
    negLog10FDR = -log10(
      pmax(FDR, .Machine$double.xmin)
    ),
    
    Significant = (
      !is.na(FDR) &
        FDR < fdr_cutoff
    ),
    
    Direction = case_when(
      log2FC > 0 ~ "Upregulated",
      log2FC < 0 ~ "Downregulated",
      TRUE ~ "Zero"
    )
  )

deg <- deg %>%
  left_join(
    gene_reference,
    by = "gene"
  ) %>%
  mutate(
    Category = ifelse(
      is.na(Category),
      "Other Genes",
      Category
    )
  )

cat("\nTotal usable genes:", nrow(deg), "\n")

cat(
  "Significant genes:",
  sum(deg$Significant),
  "\n"
)


# ============================================================
# 8. SELECT SIGNIFICANT TARGET GENES
# ============================================================

# Only genes with FDR < 0.05 are eligible.
# Positive and negative log2FC are handled separately.

eligible_targets <- deg %>%
  filter(
    Category != "Other Genes",
    Significant,
    Direction %in% c(
      "Upregulated",
      "Downregulated"
    ),
    abs(log2FC) > min_abs_log2fc,
    is.finite(negLog10FDR)
  )


# ============================================================
# 9. SELECT TOP 5 PER GROUP AND DIRECTION
# ============================================================

# Each category is processed separately for:
#   - Upregulated genes
#   - Downregulated genes
#
# Ranking:
#   1. Smallest FDR
#   2. Largest absolute log2FC as a tie-breaker

highlight_genes <- eligible_targets %>%
  group_by(
    Category,
    Direction
  ) %>%
  arrange(
    FDR,
    desc(abs(log2FC)),
    gene,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = max_genes_per_direction
  ) %>%
  ungroup()


# ============================================================
# 10. VERIFY SELECTION
# ============================================================

# All highlighted genes must pass FDR threshold.
stopifnot(
  all(highlight_genes$FDR < fdr_cutoff)
)

# No group/direction combination can exceed 5 genes.
selection_counts <- highlight_genes %>%
  count(
    Category,
    Direction,
    name = "Genes"
  )

stopifnot(
  all(
    selection_counts$Genes <= max_genes_per_direction
  )
)

cat("\n========================================\n")
cat("HIGHLIGHTED GENE COUNTS\n")
cat("========================================\n")

print(as.data.frame(selection_counts))

cat(
  "\nTotal highlighted genes:",
  nrow(highlight_genes),
  "\n"
)


# ============================================================
# 11. PREPARE BACKGROUND GENES
# ============================================================

# Every gene not selected for highlighting
# remains in the background, even if significant.

background_data <- deg %>%
  filter(
    is.finite(negLog10FDR),
    !gene %in% highlight_genes$gene
  )

cat(
  "Background genes:",
  nrow(background_data),
  "\n"
)

# ============================================================
# 12. COLORS AND SHAPES
# ============================================================

direction_colors <- c(
  "Upregulated" = "#D73027",
  "Downregulated" = "#4575B4"
)

category_shapes <- c(
  "AhR Signaling" = 16,
  "Cytokines" = 17,
  "Activation Markers" = 15
)

# ============================================================
# 13A. AXIS LIMITS
# ============================================================

# Set y-axis upper limit based on highlighted genes,
# with a little extra room for labels.
y_max_plot <- max(highlight_genes$negLog10FDR, na.rm = TRUE) + 5

# Prevent the plot from becoming too tall if one selected gene is extreme
y_max_plot <- min(y_max_plot, 40)

# Set x-axis limits based on highlighted genes plus padding
x_min_plot <- min(highlight_genes$log2FC, na.rm = TRUE) - 1
x_max_plot <- max(highlight_genes$log2FC, na.rm = TRUE) + 1

# Make the x-axis symmetric around 0
x_abs_max <- max(abs(c(x_min_plot, x_max_plot)))
x_min_plot <- -x_abs_max
x_max_plot <-  x_abs_max

# ============================================================
# 13. GENERATE TARGETED VOLCANO PLOT
# ============================================================

p <- ggplot() +
  
  # ----------------------------------------------------------
# A. BACKGROUND GENES
# ----------------------------------------------------------

geom_point(
  data = background_data,
  aes(
    x = log2FC,
    y = negLog10FDR
  ),
  color = "grey80",
  alpha = 0.45,
  size = 1.2
) +
  
  # ----------------------------------------------------------
# B. FDR REFERENCE LINE
# ----------------------------------------------------------

geom_hline(
  yintercept = -log10(fdr_cutoff),
  linetype = "dashed",
  linewidth = 0.5,
  color = "grey40"
) +
  
  # ----------------------------------------------------------
# C. ZERO FOLD CHANGE LINE
# ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  linetype = "dashed",
  linewidth = 0.5,
  color = "grey40"
) +
  
  # ----------------------------------------------------------
# D. HIGHLIGHT SELECTED GENES
# Color = direction
# Shape = biological category
# ----------------------------------------------------------

geom_point(
  data = highlight_genes,
  aes(
    x = log2FC,
    y = negLog10FDR,
    color = Direction,
    shape = Category
  ),
  size = 3.5,
  alpha = 0.95
) +
  
  # ----------------------------------------------------------
# E. GENE LABELS
# Labels colored by direction
# No connecting lines
# ----------------------------------------------------------

ggrepel::geom_text_repel(
  data = highlight_genes,
  aes(
    x = log2FC,
    y = negLog10FDR,
    label = gene,
    color = Direction
  ),
  size = 3.5,
  box.padding = 0.5,
  point.padding = 0.25,
  max.overlaps = Inf,
  max.time = 10,
  max.iter = 30000,
  seed = 123,
  segment.color = NA,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# F. DIRECTION COLORS
# ----------------------------------------------------------

scale_color_manual(
  values = direction_colors,
  name = "Regulation"
) +
  
  # ----------------------------------------------------------
# G. CATEGORY SHAPES
# ----------------------------------------------------------

scale_shape_manual(
  values = category_shapes,
  name = "Gene Category"
) +
  
  # ----------------------------------------------------------
# H. AXIS LIMITS
# ----------------------------------------------------------

coord_cartesian(
  xlim = c(x_min_plot, x_max_plot),
  ylim = c(0, y_max_plot)
) +
  
  # ----------------------------------------------------------
# I. THEME
# ----------------------------------------------------------

theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    
    plot.subtitle = element_text(
      size = 12
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    )
  ) +
  
  labs(
    title = "Indole + LPS vs LPS",
    
    subtitle = paste0(
      "AhR signaling, cytokines, and immune activation | ",
      "FDR < 0.05"
    ),
    
    x = "Log2 Fold Change (Indole + LPS / LPS)",
    
    y = expression(
      -log[10]("Adjusted P-value")
    )
  )



# ============================================================
# 14. SAVE VOLCANO PLOT
# ============================================================

# Save directly without displaying in RStudio.
# This avoids potential plotting pane viewport warnings.

ggsave(
  filename = output_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  bg = "white",
  limitsize = FALSE
)

cat("\n========================================\n")
cat("VOLCANO PLOT SAVED\n")
cat("========================================\n")

cat(output_file, "\n")


# ============================================================
# 15. PRINT HIGHLIGHTED GENES
# ============================================================

highlighted_table <- highlight_genes %>%
  select(
    gene,
    Category,
    Direction,
    log2FC,
    FDR
  ) %>%
  arrange(
    Category,
    Direction,
    FDR
  )

cat("\n========================================\n")
cat("HIGHLIGHTED GENES\n")
cat("========================================\n")

# Use a standard data.frame print without n = Inf.
print(
  as.data.frame(highlighted_table),
  row.names = FALSE
)


# ============================================================
# 16. FINISH
# ============================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "Total highlighted genes:",
  nrow(highlight_genes),
  "\n"
)

cat(
  "Maximum allowed:",
  3 * 2 * max_genes_per_direction,
  "\n"
)

cat(
  "FDR threshold:",
  fdr_cutoff,
  "\n"
)

cat(
  "Saved figure:",
  output_file,
  "\n"
)


# Gene list table
# ============================================================
# 17. SAVE GENE LISTS BY CATEGORY
# ============================================================

cat("\n========================================\n")
cat("SAVING CATEGORY GENE LISTS\n")
cat("========================================\n")

# Include all predefined genes in each category,
# including those not found in the DEG results.

category_gene_table <- gene_reference %>%
  left_join(
    deg %>%
      select(
        gene,
        log2FC,
        FDR,
        Significant,
        Direction
      ),
    by = "gene"
  ) %>%
  mutate(
    Highlighted = gene %in% highlight_genes$gene,
    Significant = coalesce(Significant, FALSE),
    Direction = ifelse(
      is.na(Direction),
      "Not Found",
      Direction
    ),
    In_DEG_Results = !is.na(log2FC)
  ) %>%
  select(
    Category,
    gene,
    In_DEG_Results,
    log2FC,
    FDR,
    Significant,
    Direction,
    Highlighted
  ) %>%
  arrange(
    Category,
    desc(Highlighted),
    FDR,
    gene
  )

# Save complete gene list
csv_file <- file.path(
  output_dir,
  "AhR_Cytokine_Activation_Gene_Lists_stats.csv"
)

write.csv(
  category_gene_table,
  csv_file,
  row.names = FALSE
)

cat("\nGene list saved to:\n", csv_file, "\n")

# Print category summary
category_summary <- category_gene_table %>%
  group_by(Category) %>%
  summarise(
    Total_Genes = n(),
    Genes_Found = sum(In_DEG_Results),
    Significant_Genes = sum(Significant),
    Highlighted_Genes = sum(Highlighted),
    .groups = "drop"
  )

print(as.data.frame(category_summary))




# ============================================================
# 17. SAVE GENE LISTS BY CATEGORY
# ============================================================

cat("\n========================================\n")
cat("SAVING CATEGORY GENE LISTS\n")
cat("========================================\n")

# Include ALL predefined genes in each category,
# regardless of expression or statistical significance.

category_gene_table <- gene_reference %>%
  group_by(Category) %>%
  summarise(
    Genes = paste(
      sort(unique(gene)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(Category)

# Save CSV
csv_file <- file.path(
  output_dir,
  "AhR_Cytokine_Activation_Gene_Lists.csv"
)

write.csv(
  category_gene_table,
  csv_file,
  row.names = FALSE
)

cat("\nGene lists saved to:\n", csv_file, "\n")

print(
  as.data.frame(category_gene_table),
  row.names = FALSE
)