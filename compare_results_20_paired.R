
#!/usr/bin/env Rscript

# ============================================================
# SALMON / DESEQ2 vs EXCEL DEG COMPARISON
#
# Excel workbook:
#   GDC_RNAseq_fold changes.xlsx
#
# Excel columns:
#   Column 1 = Gene symbol (unnamed)
#   Column 2 = Log2 fold change
#   Column 3 = -Log10 FDR p-value
#
# Outputs for EACH usable Excel sheet:
#   1. Log2FC correlation
#   2. Four-quadrant plot
#   3. DEG significance overlap
#   4. Top 20 shared DEG fold changes
#   5. Complete matched comparison CSV
#   6. Comparison statistics
#
# Also saves:
#   Master summary across all sheets
#
# FDR cutoff = 0.05
# ============================================================


# ============================================================
# 0. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "tidyr"
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

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)


# ============================================================
# 1. DIRECTORIES AND SETTINGS
# ============================================================

base_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole"

#salmon_file <- file.path(
 # base_dir,
  #"DESeq2_Paired_All20",
  #"DEG_INLPS_vs_LPS_Gene_Names.csv"
#)
salmon_file <-  "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole/DESeq2_Paired_All20/DEG_INLPS_vs_LPS_All.csv"
head(read.csv(salmon_file))
excel_file <- file.path(
  base_dir,
  "GDC_RNAseq_fold changes.xlsx"
)

output_dir <- file.path(
  base_dir,
  "Salmon_vs_Excel_Comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# ANALYSIS SETTINGS
# ------------------------------------------------------------

fdr_cutoff <- 0.05

top_n <- 20

# FALSE = Excel is INLPS vs LPS
# TRUE  = Excel is LPS vs INLPS
reverse_excel_logfc <- FALSE

# Label top discordant genes on correlation plot
label_discordant_genes <- TRUE

n_discordant_labels <- 10

# ------------------------------------------------------------
# CHECK INPUT FILES
# ------------------------------------------------------------

if (!file.exists(salmon_file)) {
  stop(
    "Salmon DEG file not found:\n",
    salmon_file
  )
}

if (!file.exists(excel_file)) {
  stop(
    "Excel workbook not found:\n",
    excel_file
  )
}

cat("\n========================================\n")
cat("SALMON VS EXCEL DEG COMPARISON\n")
cat("========================================\n")

cat("Salmon:", salmon_file, "\n")
cat("Excel:", excel_file, "\n")
cat("Output:", output_dir, "\n")


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# STANDARDIZE GENE SYMBOLS
# ------------------------------------------------------------

clean_gene <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- toupper(x)
  
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  
  x
}


# ------------------------------------------------------------
# CONVERT TO NUMERIC
# ------------------------------------------------------------

to_number <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- gsub(",", "", x, fixed = TRUE)
  
  suppressWarnings(as.numeric(x))
}


# ------------------------------------------------------------
# CONVERT -LOG10(FDR) TO FDR
# ------------------------------------------------------------

neglog_to_fdr <- function(x) {
  
  result <- 10^(-x)
  
  result[is.na(x)] <- NA_real_
  
  result
}


# ------------------------------------------------------------
# SAVE PDF AND PNG
# ------------------------------------------------------------

save_plot <- function(
    plot,
    prefix,
    width = 10,
    height = 8
) {
  
  ggsave(
    filename = paste0(prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = paste0(prefix, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}


# ------------------------------------------------------------
# SAFE DIRECTORY NAMES
# ------------------------------------------------------------

safe_filename <- function(x) {
  
  x <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    x
  )
  
  x <- gsub("^_+|_+$", "", x)
  
  x
}

# ============================================================
# 3. READ SALMON DEG RESULTS
#    CONVERT ENSEMBL IDs TO GENE SYMBOLS
# ============================================================

cat("\nReading Salmon DEG results...\n")

# ------------------------------------------------------------
# INSTALL BIOCONDUCTOR PACKAGES IF NEEDED
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

# ------------------------------------------------------------
# READ DEG RESULTS
# ------------------------------------------------------------

salmon_raw <- read.csv(
  salmon_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Original DEG rows:", nrow(salmon_raw), "\n")

# Validate required columns
required_columns <- c(
  "gene_id",
  "log2FoldChange",
  "padj"
)

missing_columns <- setdiff(
  required_columns,
  colnames(salmon_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# REMOVE ENSEMBL VERSION NUMBERS
# ------------------------------------------------------------

salmon_raw$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  salmon_raw$gene_id
)

# Example:
# ENSG00000082898.19 -> ENSG00000082898

# ------------------------------------------------------------
# MAP ENSEMBL IDs TO GENE SYMBOLS
# ------------------------------------------------------------

gene_mapping <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(salmon_raw$ensembl_id),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

salmon_raw$gene_name <- unname(
  gene_mapping[salmon_raw$ensembl_id]
)

# ------------------------------------------------------------
# REPORT MAPPING SUCCESS
# ------------------------------------------------------------

total_genes <- nrow(salmon_raw)

mapped_genes <- sum(
  !is.na(salmon_raw$gene_name) &
    salmon_raw$gene_name != ""
)

unmapped_genes <- total_genes - mapped_genes

cat("\n========================================\n")
cat("ENSEMBL GENE MAPPING RESULTS\n")
cat("========================================\n")

cat("Total genes:", total_genes, "\n")
cat("Mapped genes:", mapped_genes, "\n")
cat("Unmapped genes:", unmapped_genes, "\n")

cat(
  "Mapping percentage:",
  round(100 * mapped_genes / total_genes, 2),
  "%\n"
)

# ------------------------------------------------------------
# SAVE ANNOTATED DEG RESULTS
# ------------------------------------------------------------

annotated_file <- file.path(
  dirname(salmon_file),
  "DEG_INLPS_vs_LPS_Gene_Names.csv"
)

write.csv(
  salmon_raw,
  annotated_file,
  row.names = FALSE
)

cat(
  "\nAnnotated DEG file saved to:",
  annotated_file,
  "\n"
)

# ------------------------------------------------------------
# PREPARE DATA FOR EXCEL COMPARISON
# ------------------------------------------------------------

salmon <- salmon_raw %>%
  transmute(
    gene = clean_gene(gene_name),
    
    salmon_gene_id = as.character(gene_id),
    
    salmon_log2FC = as.numeric(log2FoldChange),
    
    salmon_FDR = as.numeric(padj),
    
    salmon_baseMean = if ("baseMean" %in% names(salmon_raw)) {
      as.numeric(baseMean)
    } else {
      NA_real_
    }
  ) %>%
  filter(
    !is.na(gene),
    is.finite(salmon_log2FC)
  ) %>%
  arrange(
    is.na(salmon_FDR),
    salmon_FDR,
    desc(salmon_baseMean)
  )

# ------------------------------------------------------------
# HANDLE DUPLICATE SYMBOLS
# ------------------------------------------------------------

duplicate_salmon <- salmon %>%
  count(gene) %>%
  filter(n > 1)

if (nrow(duplicate_salmon) > 0) {
  warning(
    nrow(duplicate_salmon),
    " duplicated gene symbols found. ",
    "Keeping the row with the smallest available FDR."
  )
}

salmon <- salmon %>%
  distinct(gene, .keep_all = TRUE)

cat("\n========================================\n")
cat("SALMON DEG SUMMARY\n")
cat("========================================\n")

cat("Unique mapped genes:", nrow(salmon), "\n")

cat(
  "Significant genes (FDR < 0.05):",
  sum(salmon$salmon_FDR < fdr_cutoff, na.rm = TRUE),
  "\n"
)
# ============================================================
# 4. READ EXCEL DEG SHEET
# ============================================================

read_excel_deg <- function(file, sheet) {
  
  cat("\nReading Excel worksheet:", sheet, "\n")
  
  dat <- readxl::read_excel(
    path = file,
    sheet = sheet,
    col_names = TRUE,
    .name_repair = "unique"
  )
  
  dat <- as.data.frame(dat)
  
  if (nrow(dat) == 0 || ncol(dat) < 3) {
    
    message(
      "Skipping sheet: fewer than 3 columns or no rows."
    )
    
    return(NULL)
  }
  
  # ==========================================================
  # IMPORTANT:
  # USE COLUMN POSITIONS INSTEAD OF COLUMN NAMES
  # ==========================================================
  
  # Column 1: Gene symbols
  gene_values <- clean_gene(dat[[1]])
  
  # Column 2: Log2 fold changes
  excel_log2FC <- to_number(dat[[2]])
  
  # Column 3: -Log10(FDR)
  excel_negLog10FDR <- to_number(dat[[3]])
  
  # ==========================================================
  # VALIDATE THE EXCEL VALUES
  # ==========================================================
  
  # A negative -log10(FDR) is not valid for FDR <= 1.
  invalid_fdr <- (
    !is.na(excel_negLog10FDR) &
      excel_negLog10FDR < 0
  )
  
  if (any(invalid_fdr)) {
    
    warning(
      "Sheet ", sheet, ": ",
      sum(invalid_fdr),
      " rows have negative -log10(FDR)."
    )
  }
  
  # Convert -log10(FDR) to FDR
  excel_FDR <- neglog_to_fdr(
    excel_negLog10FDR
  )
  
  # Reverse Excel fold changes if needed
  if (reverse_excel_logfc) {
    
    excel_log2FC <- -excel_log2FC
    
  }
  
  result <- data.frame(
    gene = gene_values,
    excel_log2FC = excel_log2FC,
    excel_negLog10FDR = excel_negLog10FDR,
    excel_FDR = excel_FDR,
    stringsAsFactors = FALSE
  )
  
  # ==========================================================
  # REMOVE INVALID ROWS
  # ==========================================================
  
  total_rows <- nrow(result)
  
  result <- result %>%
    filter(
      !is.na(gene),
      is.finite(excel_log2FC),
      !is.na(excel_negLog10FDR),
      excel_negLog10FDR >= 0,
      excel_FDR >= 0,
      excel_FDR <= 1
    )
  
  cat("Original Excel rows:", total_rows, "\n")
  cat("Valid DEG rows:", nrow(result), "\n")
  
  if (nrow(result) == 0) {
    
    message(
      "No valid DEG rows in worksheet: ",
      sheet
    )
    
    return(NULL)
  }
  
  # ==========================================================
  # HANDLE DUPLICATE GENE SYMBOLS
  # ==========================================================
  
  duplicate_excel <- result %>%
    count(gene) %>%
    filter(n > 1)
  
  if (nrow(duplicate_excel) > 0) {
    
    warning(
      "Worksheet ", sheet, " contains ",
      nrow(duplicate_excel),
      " duplicated gene symbols."
    )
    
  }
  
  # Prefer strongest evidence for duplicate symbols.
  # This avoids using FDR values that underflow to 0
  # as the ranking variable.
  result <- result %>%
    arrange(
      desc(excel_negLog10FDR),
      desc(abs(excel_log2FC))
    ) %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("Unique Excel genes:", nrow(result), "\n")
  
  return(result)
}


# ============================================================
# 5. COMPARE EACH WORKSHEET
# ============================================================

compare_sheet <- function(
    excel,
    sheet_name,
    outdir
) {
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ==========================================================
  # A. MATCH GENES
  # ==========================================================
  
  combined <- inner_join(
    salmon,
    excel,
    by = "gene"
  )
  
  if (nrow(combined) == 0) {
    
    message(
      "No matching genes in worksheet: ",
      sheet_name
    )
    
    return(NULL)
  }
  
  cat(
    "Matched genes:",
    nrow(combined),
    "\n"
  )
  
  # ==========================================================
  # B. CLASSIFY SIGNIFICANCE AND DIRECTION
  # ==========================================================
  
  combined <- combined %>%
    mutate(
      
      salmon_sig = (
        !is.na(salmon_FDR) &
          salmon_FDR < fdr_cutoff
      ),
      
      # Use the original -log10(FDR) for Excel.
      # This avoids numerical underflow.
      excel_sig = (
        excel_negLog10FDR >
          -log10(fdr_cutoff)
      ),
      
      salmon_direction = case_when(
        salmon_log2FC > 0 ~ "Up",
        salmon_log2FC < 0 ~ "Down",
        TRUE ~ "Zero"
      ),
      
      excel_direction = case_when(
        excel_log2FC > 0 ~ "Up",
        excel_log2FC < 0 ~ "Down",
        TRUE ~ "Zero"
      ),
      
      direction_agreement = case_when(
        salmon_log2FC * excel_log2FC > 0 ~ "Agree",
        salmon_log2FC * excel_log2FC < 0 ~ "Disagree",
        TRUE ~ "Zero fold change"
      ),
      
      significance_group = case_when(
        salmon_sig & excel_sig ~ "Both Significant",
        salmon_sig & !excel_sig ~ "Salmon Only",
        !salmon_sig & excel_sig ~ "Excel Only",
        TRUE ~ "Neither"
      ),
      
      both_significant = salmon_sig & excel_sig,
      
      quadrant = case_when(
        excel_log2FC > 0 &
          salmon_log2FC > 0 ~ "Up in Both",
        
        excel_log2FC < 0 &
          salmon_log2FC < 0 ~ "Down in Both",
        
        excel_log2FC > 0 &
          salmon_log2FC < 0 ~ "Excel Up / Salmon Down",
        
        excel_log2FC < 0 &
          salmon_log2FC > 0 ~ "Excel Down / Salmon Up",
        
        TRUE ~ "Zero"
      ),
      
      absolute_FC_difference = abs(
        salmon_log2FC - excel_log2FC
      )
    )
  
  # ==========================================================
  # C. SAVE COMPLETE COMPARISON TABLE
  # ==========================================================
  
  combined <- combined %>%
    arrange(
      desc(both_significant),
      salmon_FDR
    )
  
  write.csv(
    combined,
    file.path(
      outdir,
      "Combined_Salmon_Excel_DEGs.csv"
    ),
    row.names = FALSE
  )
  
  # Also save genes with opposite directions
  discordant <- combined %>%
    filter(
      direction_agreement == "Disagree"
    ) %>%
    arrange(desc(absolute_FC_difference))
  
  write.csv(
    discordant,
    file.path(
      outdir,
      "Discordant_Genes.csv"
    ),
    row.names = FALSE
  )
  
  # ==========================================================
  # D. CORRELATION STATISTICS
  # ==========================================================
  
  cor_dat <- combined %>%
    filter(
      is.finite(salmon_log2FC),
      is.finite(excel_log2FC)
    )
  
  pearson_r <- if (
    nrow(cor_dat) >= 3 &&
    sd(cor_dat$salmon_log2FC) > 0 &&
    sd(cor_dat$excel_log2FC) > 0
  ) {
    
    cor(
      cor_dat$excel_log2FC,
      cor_dat$salmon_log2FC,
      method = "pearson"
    )
    
  } else {
    
    NA_real_
    
  }
  
  spearman_r <- if (
    nrow(cor_dat) >= 3 &&
    sd(cor_dat$salmon_log2FC) > 0 &&
    sd(cor_dat$excel_log2FC) > 0
  ) {
    
    suppressWarnings(
      cor(
        cor_dat$excel_log2FC,
        cor_dat$salmon_log2FC,
        method = "spearman"
      )
    )
    
  } else {
    
    NA_real_
    
  }
  
  # ==========================================================
  # E. LOG2FC CORRELATION PLOT
  # ==========================================================
  
  p_cor <- ggplot(
    cor_dat,
    aes(
      x = excel_log2FC,
      y = salmon_log2FC
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      color = "black"
    ) +
    
    geom_point(
      aes(color = direction_agreement),
      alpha = 0.6,
      size = 1.8
    ) +
    
    scale_color_manual(
      values = c(
        "Agree" = "#3182BD",
        "Disagree" = "#E6550D",
        "Zero fold change" = "grey70"
      )
    ) +
    
    coord_equal() +
    
    theme_classic(base_size = 13) +
    
    labs(
      title = paste(
        "Log2FC Correlation:",
        sheet_name
      ),
      
      subtitle = sprintf(
        "Pearson r = %.3f | Spearman rho = %.3f | n = %d",
        pearson_r,
        spearman_r,
        nrow(cor_dat)
      ),
      
      x = "Excel Log2FC",
      y = "Salmon Log2FC",
      
      color = "Direction Agreement"
    )
  
  # Label discordant genes with largest FC differences
  if (
    label_discordant_genes &&
    nrow(discordant) > 0
  ) {
    
    label_dat <- discordant %>%
      slice_head(
        n = n_discordant_labels
      )
    
    p_cor <- p_cor +
      
      geom_text_repel(
        data = label_dat,
        aes(label = gene),
        size = 3,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      )
    
  }
  
  save_plot(
    p_cor,
    file.path(
      outdir,
      "01_Log2FC_Correlation"
    ),
    width = 11,
    height = 9
  )
  
  # ==========================================================
  # F. FOUR-QUADRANT DEG PLOT
  # ==========================================================
  
  # Keep only genes significant in both analyses (FDR < 0.05)
  quadrant_data <- combined %>%
    filter(both_significant)
  
  quadrant_counts <- quadrant_data %>%
    count(
      quadrant,
      name = "Genes"
    )
  
  write.csv(
    quadrant_counts,
    file.path(
      outdir,
      "Four_Quadrant_Counts.csv"
    ),
    row.names = FALSE
  )
  
  p_quad <- ggplot(
    quadrant_data,
    aes(
      x = excel_log2FC,
      y = salmon_log2FC
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      color = "grey40"
    ) +
    
    geom_vline(
      xintercept = 0,
      color = "grey40"
    ) +
    
    geom_point(
      aes(color = quadrant),
      alpha = 0.65,
      size = 1.7
    ) +
    
    scale_color_manual(
      values = c(
        "Up in Both" = "#D73027",
        "Down in Both" = "#4575B4",
        "Excel Up / Salmon Down" = "#FDAE61",
        "Excel Down / Salmon Up" = "#91BFDB",
        "Zero" = "grey70"
      )
    ) +
    
    theme_classic(base_size = 13) +
    
    labs(
      title = paste(
        "Four-Quadrant DEG Comparison:",
        sheet_name
      ),
      
      subtitle = paste0(
        nrow(quadrant_data),
        " genes significant in both analyses (FDR < 0.05)"
      ),
      
      x = "Excel Log2FC",
      y = "Salmon Log2FC",
      
      color = "Direction"
    )
  
  save_plot(
    p_quad,
    file.path(
      outdir,
      "02_Four_Quadrant_DEGs"
    ),
    width = 11,
    height = 9
  )
  
  # ==========================================================
  # G. SIGNIFICANT DEG OVERLAP
  # ==========================================================
  
  overlap_levels <- c(
    "Both Significant",
    "Salmon Only",
    "Excel Only",
    "Neither"
  )
  
  overlap_counts <- combined %>%
    count(
      significance_group,
      name = "Genes"
    ) %>%
    mutate(
      significance_group = factor(
        significance_group,
        levels = overlap_levels
      )
    ) %>%
    tidyr::complete(
      significance_group,
      fill = list(Genes = 0)
    )
  
  write.csv(
    overlap_counts,
    file.path(
      outdir,
      "DEG_Overlap_Counts.csv"
    ),
    row.names = FALSE
  )
  
  p_overlap <- ggplot(
    overlap_counts,
    aes(
      x = significance_group,
      y = Genes,
      fill = significance_group
    )
  ) +
    
    geom_col(width = 0.65) +
    
    geom_text(
      aes(label = Genes),
      vjust = -0.4,
      size = 4
    ) +
    
    scale_fill_manual(
      values = c(
        "Both Significant" = "#238B45",
        "Salmon Only" = "#3182BD",
        "Excel Only" = "#E6550D",
        "Neither" = "grey65"
      )
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.12)
      )
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 15,
        hjust = 1
      )
    ) +
    
    labs(
      title = paste(
        "DEG Overlap:",
        sheet_name
      ),
      
      subtitle = paste0(
        "FDR < 0.05 | Among ",
        nrow(combined),
        " matched genes"
      ),
      
      x = NULL,
      y = "Number of Genes"
    )
  
  save_plot(
    p_overlap,
    file.path(
      outdir,
      "03_DEG_Overlap"
    ),
    width = 10,
    height = 8
  )
  
  # ==========================================================
  # H. TOP 20 SHARED SIGNIFICANT DEGS
  # ==========================================================
  
  shared_sig <- combined %>%
    filter(both_significant)
  
  # Rank by statistical evidence in both analyses.
  # Use Excel's original -log10(FDR) to avoid underflow.
  
  top_shared <- shared_sig %>%
    mutate(
      combined_score =
        -log10(
          pmax(
            salmon_FDR,
            .Machine$double.xmin
          )
        ) +
        excel_negLog10FDR
    ) %>%
    arrange(desc(combined_score)) %>%
    slice_head(n = top_n)
  
  write.csv(
    top_shared,
    file.path(
      outdir,
      "Top_Shared_Significant_Genes.csv"
    ),
    row.names = FALSE
  )
  
  if (nrow(top_shared) > 0) {
    
    top_long <- top_shared %>%
      select(
        gene,
        salmon_log2FC,
        excel_log2FC
      ) %>%
      pivot_longer(
        cols = c(
          salmon_log2FC,
          excel_log2FC
        ),
        
        names_to = "Analysis",
        values_to = "Log2FC"
      ) %>%
      mutate(
        Analysis = recode(
          Analysis,
          salmon_log2FC = "Salmon",
          excel_log2FC = "Excel"
        ),
        
        gene = factor(
          gene,
          levels = rev(top_shared$gene)
        )
      )
    
    p_top <- ggplot(
      top_long,
      aes(
        x = gene,
        y = Log2FC,
        fill = Analysis
      )
    ) +
      
      geom_col(
        position = position_dodge(
          width = 0.75
        ),
        width = 0.7
      ) +
      
      geom_hline(
        yintercept = 0,
        color = "grey40"
      ) +
      
      coord_flip() +
      
      scale_fill_manual(
        values = c(
          "Salmon" = "#3182BD",
          "Excel" = "#E6550D"
        )
      ) +
      
      theme_classic(base_size = 13) +
      
      labs(
        title = paste(
          "Top Shared DEGs:",
          sheet_name
        ),
        
        subtitle = paste0(
          "Top ",
          nrow(top_shared),
          " genes significant in both analyses"
        ),
        
        x = NULL,
        y = "Log2 Fold Change",
        fill = "Analysis"
      )
    
    save_plot(
      p_top,
      file.path(
        outdir,
        "04_Top_Shared_DEG_Log2FC"
      ),
      width = 11,
      height = 9
    )
    
  } else {
    
    message(
      "No genes significant in both analyses: ",
      sheet_name
    )
    
  }
  
  # ==========================================================
  # I. SUMMARY STATISTICS
  # ==========================================================
  
  agreement_n <- sum(
    combined$direction_agreement == "Agree"
  )
  
  shared_agreement <- if (
    nrow(shared_sig) > 0
  ) {
    
    100 * mean(
      shared_sig$direction_agreement == "Agree"
    )
    
  } else {
    
    NA_real_
    
  }
  
  summary <- data.frame(
    
    worksheet = sheet_name,
    
    excel_genes = nrow(excel),
    
    salmon_genes = nrow(salmon),
    
    matched_genes = nrow(combined),
    
    pearson_r = pearson_r,
    
    spearman_rho = spearman_r,
    
    direction_agreement_pct =
      100 * agreement_n / nrow(combined),
    
    both_significant = nrow(shared_sig),
    
    salmon_only = sum(
      combined$significance_group == "Salmon Only"
    ),
    
    excel_only = sum(
      combined$significance_group == "Excel Only"
    ),
    
    neither = sum(
      combined$significance_group == "Neither"
    ),
    
    both_sig_direction_agreement_pct =
      shared_agreement,
    
    discordant_genes = nrow(discordant)
    
  )
  
  write.csv(
    summary,
    file.path(
      outdir,
      "Comparison_Summary.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nWorksheet summary:\n")
  print(summary)
  
  return(summary)
}


# ============================================================
# 6. PROCESS EVERY EXCEL WORKSHEET
# ============================================================

sheets <- readxl::excel_sheets(excel_file)

cat("\nExcel worksheets found:\n")
print(sheets)

all_summaries <- list()

for (i in seq_along(sheets)) {
  
  sheet <- sheets[i]
  
  cat("\n========================================\n")
  cat("PROCESSING:", sheet, "\n")
  cat("========================================\n")
  
  excel <- read_excel_deg(
    file = excel_file,
    sheet = sheet
  )
  
  if (is.null(excel)) {
    next
  }
  
  safe_name <- paste0(
    sprintf("%02d", i),
    "_",
    safe_filename(sheet)
  )
  
  sheet_dir <- file.path(
    output_dir,
    safe_name
  )
  
  summary <- compare_sheet(
    excel = excel,
    sheet_name = sheet,
    outdir = sheet_dir
  )
  
  if (!is.null(summary)) {
    
    all_summaries[[length(all_summaries) + 1]] <- summary
    
  }
}


# ============================================================
# 7. MASTER SUMMARY
# ============================================================

if (length(all_summaries) == 0) {
  
  stop(
    "No usable DEG worksheets found. ",
    "Confirm columns 1-3 contain gene symbols, ",
    "log2FC, and -log10(FDR)."
  )
  
}

master_summary <- bind_rows(
  all_summaries
)

write.csv(
  master_summary,
  file.path(
    output_dir,
    "All_Worksheet_Comparison_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("ALL COMPARISONS COMPLETE\n")
cat("========================================\n")

cat(
  "Worksheets analyzed:",
  nrow(master_summary),
  "\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)

print(master_summary)

cat("\nDONE.\n")




# ============================================================
# 9. ADDITIONAL COMPARISON STATISTICS
# ============================================================
#
# Produces:
#   1. Overall_Comparison_Statistics.csv
#   2. Gene_Level_Comparison_Table.csv
#
# Uses all processed Excel worksheets.
#
# Similar log2FC:
#   Absolute difference <= 0.5
#
# Same direction:
#   Both log2FC values positive or both negative
# ============================================================

cat("\n========================================\n")
cat("ADDITIONAL COMPARISON STATISTICS\n")
cat("========================================\n")

fc_similarity_threshold <- 0.5

overall_results <- list()
gene_level_results <- list()

for (i in seq_along(sheets)) {
  
  sheet <- sheets[i]
  
  cat("\nSummarizing worksheet:", sheet, "\n")
  
  # Read the Excel worksheet
  excel <- read_excel_deg(
    file = excel_file,
    sheet = sheet
  )
  
  if (is.null(excel)) next
  
  # ==========================================================
  # A. IDENTIFY MATCHING GENES
  # ==========================================================
  
  combined <- inner_join(
    salmon,
    excel,
    by = "gene"
  )
  
  if (nrow(combined) == 0) next
  
  # ==========================================================
  # B. CLASSIFY EACH GENE
  # ==========================================================
  
  combined <- combined %>%
    mutate(
      
      salmon_significant = (
        !is.na(salmon_FDR) &
          salmon_FDR < fdr_cutoff
      ),
      
      excel_significant = (
        !is.na(excel_negLog10FDR) &
          excel_negLog10FDR > -log10(fdr_cutoff)
      ),
      
      significant_in_both = (
        salmon_significant &
          excel_significant
      ),
      
      same_direction = (
        (salmon_log2FC > 0 & excel_log2FC > 0) |
          (salmon_log2FC < 0 & excel_log2FC < 0)
      ),
      
      log2FC_difference = abs(
        salmon_log2FC - excel_log2FC
      ),
      
      similar_log2FC = (
        log2FC_difference <= fc_similarity_threshold
      ),
      
      same_direction_and_similar = (
        same_direction & similar_log2FC
      ),
      
      significant_and_concordant = (
        significant_in_both &
          same_direction &
          similar_log2FC
      ),
      
      comparison_category = case_when(
        significant_and_concordant ~
          "Significant in both, same direction, similar FC",
        
        significant_in_both & same_direction ~
          "Significant in both, same direction, different FC",
        
        significant_in_both & !same_direction ~
          "Significant in both, opposite direction",
        
        salmon_significant & !excel_significant ~
          "Significant in Salmon only",
        
        !salmon_significant & excel_significant ~
          "Significant in Excel only",
        
        TRUE ~ "Not significant in either"
      ),
      
      worksheet = sheet
    )
  
  # ==========================================================
  # C. GENE-LEVEL COMPARISON TABLE
  # ==========================================================
  
  gene_level <- combined %>%
    select(
      worksheet,
      gene,
      salmon_gene_id,
      salmon_log2FC,
      excel_log2FC,
      log2FC_difference,
      salmon_FDR,
      excel_FDR,
      excel_negLog10FDR,
      salmon_significant,
      excel_significant,
      significant_in_both,
      same_direction,
      similar_log2FC,
      same_direction_and_similar,
      significant_and_concordant,
      comparison_category
    ) %>%
    arrange(
      desc(significant_and_concordant),
      desc(significant_in_both),
      log2FC_difference
    )
  
  gene_level_results[[length(gene_level_results) + 1]] <-
    gene_level
  
  # ==========================================================
  # D. OVERALL STATISTICS
  # ==========================================================
  
  salmon_sig_total <- sum(
    !is.na(salmon$salmon_FDR) &
      salmon$salmon_FDR < fdr_cutoff
  )
  
  excel_sig_total <- sum(
    excel$excel_negLog10FDR > -log10(fdr_cutoff),
    na.rm = TRUE
  )
  
  n_matched <- nrow(combined)
  
  n_both_sig <- sum(combined$significant_in_both)
  n_same_direction <- sum(combined$same_direction)
  n_similar_fc <- sum(combined$similar_log2FC)
  
  n_concordant <- sum(
    combined$same_direction_and_similar
  )
  
  n_sig_concordant <- sum(
    combined$significant_and_concordant
  )
  
  n_both_sig_same_direction <- sum(
    combined$significant_in_both &
      combined$same_direction
  )
  
  n_both_sig_similar_fc <- sum(
    combined$significant_in_both &
      combined$similar_log2FC
  )
  
  overall <- data.frame(
    
    Worksheet = sheet,
    
    Salmon_Total_Genes = nrow(salmon),
    
    Excel_Total_Genes = nrow(excel),
    
    Matched_Genes = n_matched,
    
    Salmon_Significant_FDR005 = salmon_sig_total,
    
    Excel_Significant_FDR005 = excel_sig_total,
    
    Significant_In_Both = n_both_sig,
    
    Same_Direction = n_same_direction,
    
    Similar_Log2FC = n_similar_fc,
    
    Same_Direction_And_Similar_Log2FC = n_concordant,
    
    Significant_Both_Same_Direction =
      n_both_sig_same_direction,
    
    Significant_Both_Similar_Log2FC =
      n_both_sig_similar_fc,
    
    Significant_Both_Same_Direction_And_Similar_Log2FC =
      n_sig_concordant,
    
    Matched_Direction_Agreement_Percent =
      round(100 * n_same_direction / n_matched, 2),
    
    Matched_Overall_Concordance_Percent =
      round(100 * n_concordant / n_matched, 2),
    
    Significant_Both_Direction_Agreement_Percent =
      if (n_both_sig > 0) {
        round(
          100 * n_both_sig_same_direction / n_both_sig,
          2
        )
      } else {
        NA_real_
      },
    
    Significant_Both_FC_Concordance_Percent =
      if (n_both_sig > 0) {
        round(
          100 * n_sig_concordant / n_both_sig,
          2
        )
      } else {
        NA_real_
      },
    
    stringsAsFactors = FALSE
  )
  
  overall_results[[length(overall_results) + 1]] <-
    overall
}

# ============================================================
# E. COMBINE ALL WORKSHEETS
# ============================================================

if (length(overall_results) == 0) {
  stop("No valid worksheet comparisons were available.")
}

overall_table <- bind_rows(overall_results)

gene_level_table <- bind_rows(gene_level_results)

# ============================================================
# F. SAVE TABLES
# ============================================================

write.csv(
  overall_table,
  file.path(
    output_dir,
    "Overall_Comparison_Statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_level_table,
  file.path(
    output_dir,
    "Gene_Level_Comparison_Table.csv"
  ),
  row.names = FALSE
)

# ============================================================
# G. PRINT SUMMARY
# ============================================================

cat("\n========================================\n")
cat("OVERALL COMPARISON RESULTS\n")
cat("========================================\n")

print(
  overall_table,
  row.names = FALSE
)

cat("\nComparison threshold:\n")
cat(
  "Similar log2FC = absolute difference <=",
  fc_similarity_threshold,
  "\n"
)

cat("\nTables saved to:\n")
cat(output_dir, "\n")

cat("\nADDITIONAL COMPARISON COMPLETE\n")





# anno volcano
# ============================================================
# 3B. ANNOTATED SALMON / DESEQ2 VOLCANO PLOT
# ============================================================

cat("\nGenerating annotated volcano plot...\n")

volcano <- salmon %>%
  mutate(
    negLog10FDR = -log10(
      pmax(salmon_FDR, .Machine$double.xmin)
    ),
    
    significance = case_when(
      !is.na(salmon_FDR) &
        salmon_FDR < 0.1 &
        salmon_log2FC > 0 ~ "Upregulated",
      
      !is.na(salmon_FDR) &
        salmon_FDR < 0.1 &
        salmon_log2FC < 0 ~ "Downregulated",
      
      TRUE ~ "Not Significant"
    )
  )

# Select top 15 genes in each direction
top_up <- volcano %>%
  filter(significance == "Upregulated") %>%
  arrange(salmon_FDR, desc(salmon_log2FC)) %>%
  slice_head(n = 15)

top_down <- volcano %>%
  filter(significance == "Downregulated") %>%
  arrange(salmon_FDR, salmon_log2FC) %>%
  slice_head(n = 15)

label_genes <- bind_rows(top_up, top_down)

# Generate plot
p_volcano <- ggplot(
  volcano,
  aes(
    x = salmon_log2FC,
    y = negLog10FDR
  )
) +
  geom_point(
    aes(color = significance),
    alpha = 0.6,
    size = 1.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    seed = 123
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Not Significant" = "grey75"
    )
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Indole + LPS vs LPS",
    subtitle = "Salmon / DESeq2 | Top 15 genes per direction",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "DEG Status"
  )

# Save annotated volcano
ggsave(
  filename = file.path(
    output_dir,
    "Volcano_INLPS_vs_LPS_Annotated.pdf"
  ),
  plot = p_volcano,
  width = 12,
  height = 10
)

ggsave(
  filename = file.path(
    output_dir,
    "Volcano_INLPS_vs_LPS_Annotated.png"
  ),
  plot = p_volcano,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Annotated volcano plot saved.\n")




# anno bar
# ============================================================
# 10. BAR CHART: VOLCANO-ANNOTATED GENES
#     SALMON / DESEQ2 vs EXCEL
#
# Selects:
#   - Top 15 significant upregulated Salmon genes
#   - Top 15 significant downregulated Salmon genes
#
# Ranking:
#   - Smallest Salmon adjusted p-value (FDR)
#
# For each Excel worksheet:
#   - Matches the selected genes by gene symbol
#   - Compares Salmon and Excel log2FC
#   - Saves PDF, PNG, and CSV
#
# IMPORTANT:
# Uses the same selection criteria as the annotated
# Salmon volcano plot from Section 3B.
# ============================================================

cat("\n========================================\n")
cat("VOLCANO-ANNOTATED GENE BAR CHARTS\n")
cat("========================================\n")

# ------------------------------------------------------------
# 1. SELECT THE SAME GENES AS THE VOLCANO PLOT
# ------------------------------------------------------------

volcano_genes <- salmon %>%
  mutate(
    significance = case_when(
      !is.na(salmon_FDR) &
        salmon_FDR < fdr_cutoff &
        salmon_log2FC > 0 ~ "Upregulated",
      
      !is.na(salmon_FDR) &
        salmon_FDR < fdr_cutoff &
        salmon_log2FC < 0 ~ "Downregulated",
      
      TRUE ~ "Not Significant"
    )
  )

# Top 15 upregulated
top_up <- volcano_genes %>%
  filter(significance == "Upregulated") %>%
  arrange(salmon_FDR, desc(salmon_log2FC)) %>%
  slice_head(n = 15)

# Top 15 downregulated
top_down <- volcano_genes %>%
  filter(significance == "Downregulated") %>%
  arrange(salmon_FDR, salmon_log2FC) %>%
  slice_head(n = 15)

# Combine annotated genes
volcano_labels <- bind_rows(
  top_up,
  top_down
)

cat(
  "Selected volcano genes:",
  nrow(volcano_labels),
  "\n"
)

# Save selected gene list
write.csv(
  volcano_labels,
  file.path(
    output_dir,
    "Volcano_Annotated_Gene_List.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. COMPARE TO EACH EXCEL WORKSHEET
# ------------------------------------------------------------

for (i in seq_along(sheets)) {
  
  sheet <- sheets[i]
  
  cat("\nProcessing worksheet:", sheet, "\n")
  
  # Read Excel using existing function
  excel <- read_excel_deg(
    file = excel_file,
    sheet = sheet
  )
  
  if (is.null(excel)) next
  
  # Match selected volcano genes with Excel
  combined <- volcano_labels %>%
    inner_join(
      excel,
      by = "gene"
    )
  
  cat(
    "Matched volcano genes:",
    nrow(combined),
    "of",
    nrow(volcano_labels),
    "\n"
  )
  
  if (nrow(combined) == 0) next
  
  # ----------------------------------------------------------
  # 3. ADD GENE-LEVEL COMPARISON STATISTICS
  # ----------------------------------------------------------
  
  combined <- combined %>%
    mutate(
      salmon_significant = (
        !is.na(salmon_FDR) &
          salmon_FDR < fdr_cutoff
      ),
      
      excel_significant = (
        !is.na(excel_negLog10FDR) &
          excel_negLog10FDR > -log10(fdr_cutoff)
      ),
      
      significant_in_both = (
        salmon_significant & excel_significant
      ),
      
      same_direction = (
        (salmon_log2FC > 0 & excel_log2FC > 0) |
          (salmon_log2FC < 0 & excel_log2FC < 0)
      ),
      
      log2FC_difference = abs(
        salmon_log2FC - excel_log2FC
      ),
      
      similar_log2FC = log2FC_difference <= 0.5
    )
  
  # ----------------------------------------------------------
  # 4. SAVE COMPARISON TABLE
  # ----------------------------------------------------------
  
  sheet_dir <- file.path(
    output_dir,
    paste0(
      sprintf("%02d", i),
      "_",
      safe_filename(sheet)
    )
  )
  
  dir.create(
    sheet_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    combined,
    file.path(
      sheet_dir,
      "Volcano_Annotated_Genes_Comparison.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 5. FORMAT DATA FOR GROUPED BAR CHART
  # ----------------------------------------------------------
  
  # Sort by Salmon log2FC
  gene_order <- combined %>%
    arrange(salmon_log2FC) %>%
    pull(gene)
  
  plot_data <- combined %>%
    select(
      gene,
      salmon_log2FC,
      excel_log2FC
    ) %>%
    pivot_longer(
      cols = c(
        salmon_log2FC,
        excel_log2FC
      ),
      names_to = "Analysis",
      values_to = "Log2FC"
    ) %>%
    mutate(
      Analysis = recode(
        Analysis,
        salmon_log2FC = "Salmon",
        excel_log2FC = "Excel"
      ),
      
      gene = factor(
        gene,
        levels = gene_order
      )
    )
  
  # ----------------------------------------------------------
  # 6. GENERATE BAR CHART
  # ----------------------------------------------------------
  
  p <- ggplot(
    plot_data,
    aes(
      x = gene,
      y = Log2FC,
      fill = Analysis
    )
  ) +
    
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.7
    ) +
    
    geom_hline(
      yintercept = 0,
      color = "grey40",
      linewidth = 0.4
    ) +
    
    coord_flip() +
    
    scale_fill_manual(
      values = c(
        "Salmon" = "#3182BD",
        "Excel" = "#E6550D"
      )
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold")
    ) +
    
    labs(
      title = paste(
        "Volcano-Annotated Genes:",
        sheet
      ),
      
      subtitle = paste0(
        nrow(combined),
        " of ",
        nrow(volcano_labels),
        " volcano-annotated genes matched"
      ),
      
      x = NULL,
      y = "Log2 Fold Change",
      fill = "Analysis"
    )
  
  # ----------------------------------------------------------
  # 7. SAVE PDF AND PNG
  # ----------------------------------------------------------
  
  save_plot(
    p,
    file.path(
      sheet_dir,
      "05_Volcano_Annotated_Genes_Log2FC_Comparison"
    ),
    width = 11,
    height = 11
  )
  
  # ----------------------------------------------------------
  # 8. PRINT SUMMARY
  # ----------------------------------------------------------
  
  cat(
    "Significant in both:",
    sum(combined$significant_in_both),
    "\n"
  )
  
  cat(
    "Same direction:",
    sum(combined$same_direction),
    "\n"
  )
  
  cat(
    "Similar log2FC (difference <= 0.5):",
    sum(combined$similar_log2FC),
    "\n"
  )
  
  cat(
    "Same direction and similar log2FC:",
    sum(
      combined$same_direction &
        combined$similar_log2FC
    ),
    "\n"
  )
}

cat("\n========================================\n")
cat("VOLCANO GENE COMPARISONS COMPLETE\n")
cat("========================================\n")