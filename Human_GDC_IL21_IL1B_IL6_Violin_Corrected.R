#!/usr/bin/env Rscript
# ============================================================
# HUMAN GDC: IL21, IL1B, IL6 EXPRESSION VIOLIN PLOTS
# Indole + LPS versus LPS alone; 10 matched sample pairs
# Input: DESeq2_Paired_All20/Normalized_Gene_Counts.csv
#        unnamed first column containing versioned Ensembl IDs
# Outputs: three separate PNG/PDF figures plus audit CSVs
# ============================================================

base_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Human_GDC_Indole"
counts_file <- file.path(base_dir, "DESeq2_Paired_All20", "Normalized_Gene_Counts.csv")
output_dir <- file.path(base_dir, "Gene_Violin_Plots")
show_paired_lines <- TRUE # ONLY if consecutive samples are biologically paired

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

if (!file.exists(counts_file)) stop("Input not found: ", counts_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Use the exact sample columns and labels supplied for this study.
samples <- data.frame(
  sample = c(
    "1_GDC_indlps_11_7_2024", "2_GDC_lps_11_7_2024",
    "3_169936_inlps_11_25_2024", "4_169936_lps_11_25_2024",
    "5_001716_inlps_11_25_2024", "6_001716_lps_11_25_2024",
    "7_169968_inlps_11_25_2024", "8_169968_lps_11_25_2024",
    "9_001528_in_plus_lps_12_12_2025", "10_001528_lps_12_12_2024",
    "11_12_13_2024_No11", "12_12_13_2024_No12",
    "13_12_16_2024_No9", "14_12_16_2024_No10",
    "15_1_9_2025_No3", "16_1_9_2025_No4",
    "17_1_9_2025_No8", "18_1_9_2025_No9",
    "19_1_9_2025_No13", "20_1_9_2025_No14"
  ),
  treatment = rep(c("Indole + LPS", "LPS"), 10),
  pair_id = rep(sprintf("Pair_%02d", 1:10), each = 2),
  stringsAsFactors = FALSE
)
samples$treatment <- factor(samples$treatment, levels = c("LPS", "Indole + LPS"))

# Preserve the unnamed first column rather than assuming 'gene_id' exists.
counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts) < 2) stop("Normalized counts file has fewer than 2 columns")
if (!"gene_id" %in% names(counts)) {
  first_name <- names(counts)[1]
  first_values <- trimws(as.character(counts[[1]]))
  is_ensg <- !is.na(first_values) & grepl("^ENSG[0-9]+(\\.[0-9]+)?$", first_values)
  if (!any(is_ensg)) {
    stop("First column (named '", first_name, "') does not appear to contain Ensembl IDs. First values: ",
         paste(head(first_values, 5), collapse = ", "))
  }
  names(counts)[1] <- "gene_id"
}
counts$ensembl_id <- sub("\\.[0-9]+$", "", trimws(as.character(counts$gene_id)))

missing_samples <- setdiff(samples$sample, names(counts))
if (length(missing_samples)) {
  stop("Sample columns missing from normalized counts: ", paste(missing_samples, collapse = ", "))
}
cat("Loaded", nrow(counts), "gene rows and all", nrow(samples), "sample columns.\n")
cat("First Ensembl IDs:", paste(head(counts$ensembl_id, 4), collapse = ", "), "\n")

# Stable human Ensembl IDs, without version suffixes.
target_genes <- c(IL21 = "ENSG00000138684",
                  IL1B = "ENSG00000125538",
                  IL6  = "ENSG00000136244")

write.csv(samples, file.path(output_dir, "Sample_Treatment_Mapping.csv"), row.names = FALSE)
all_expression <- list()
lookup <- data.frame(gene = names(target_genes), ensembl_id = unname(target_genes),
                     n_rows_found = integer(length(target_genes)))

for (gene in names(target_genes)) {
  matching_rows <- which(counts$ensembl_id == target_genes[[gene]])
  lookup$n_rows_found[lookup$gene == gene] <- length(matching_rows)
  cat(gene, "(", target_genes[[gene]], "): ", length(matching_rows), " row(s)\n", sep = "")
  if (length(matching_rows) == 0) {
    warning("Gene not found in normalized count matrix: ", gene)
    next
  }
  if (length(matching_rows) != 1) {
    stop("Multiple rows match ", gene, ". Resolve duplicated Ensembl IDs before plotting.")
  }
  raw_values <- unlist(counts[matching_rows, samples$sample, drop = FALSE], use.names = FALSE)
  values <- suppressWarnings(as.numeric(raw_values))
  if (length(values) != nrow(samples) || any(!is.finite(values)) || any(values < 0)) {
    stop("Invalid normalized counts for ", gene, ": ensure 20 nonnegative numeric sample values")
  }
  plot_df <- samples
  plot_df$gene <- gene
  plot_df$ensembl_id <- target_genes[[gene]]
  plot_df$normalized_count <- values
  plot_df$log2_normalized_count <- log2(values + 1)
  all_expression[[gene]] <- plot_df

  # Violin distribution, embedded boxplot, individual sample points,
  # optional pair connections (valid only if stated pairings are biological).
  p <- ggplot(plot_df, aes(x = treatment, y = log2_normalized_count)) +
    geom_violin(aes(fill = treatment), trim = FALSE, alpha = 0.45,
                color = "grey30", width = 0.8) +
    geom_boxplot(width = 0.17, outlier.shape = NA, fill = "white",
                 alpha = 0.8, linewidth = 0.5)
  if (show_paired_lines) {
    p <- p + geom_line(aes(group = pair_id), color = "grey50",
                       linewidth = 0.45, alpha = 0.45)
  }
  p <- p +
    geom_point(aes(fill = treatment), shape = 21, color = "grey15",
               size = 2.4, stroke = 0.35, alpha = 0.9,
               position = position_jitter(width = 0.045, height = 0, seed = 123)) +
    scale_fill_manual(values = c("LPS" = "#4575B4", "Indole + LPS" = "#D73027")) +
    labs(title = paste0(gene, " expression"),
         subtitle = "Human GDC | 10 samples per treatment",
         x = NULL, y = "log2(DESeq2 normalized counts + 1)") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none", axis.text.x = element_text(size = 12))

  prefix <- file.path(output_dir, paste0(gene, "_Violin_IndoleLPS_vs_LPS"))
  ggsave(paste0(prefix, ".png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot = p, width = 7, height = 6)
  cat("Saved ", gene, " violin to ", output_dir, "\n", sep = "")
}

write.csv(lookup, file.path(output_dir, "Target_Gene_Lookup.csv"), row.names = FALSE)
if (!length(all_expression)) stop("None of the target genes were found; inspect Target_Gene_Lookup.csv")
write.csv(do.call(rbind, all_expression),
          file.path(output_dir, "IL21_IL1B_IL6_Normalized_Expression.csv"), row.names = FALSE)
cat("DONE. Results saved to: ", output_dir, "\n", sep = "")
cat("NOTE: Connections are based on consecutive samples; verify that these are real biological pairs.\n")
