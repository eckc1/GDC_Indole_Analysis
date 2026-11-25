#!/usr/bin/env Rscript
# ============================================================
# SALMON -> QC -> PAIRED DESEQ2: ALL 20 SAMPLES / 10 PAIRS
# Indole + LPS versus LPS (reference)
#
# Sample IDs and treatments: 'samples for RNA.xlsx', Sheets 2 & 5
# Pair definitions: Sheet 3 (1-2, 3-4, ..., 19-20)
# Salmon directory names: derived from Sheet 5 FASTQ basenames,
# dropping _S##_L008_R#_001.fastq(.gz) suffixes.
#
# Assumptions to VERIFY with lab:
# - Each indicated pair represents the SAME biological donor.
# - Sheet1 calls sample 7 donor 169368, but sample labels say 169968.
# - Sample 9 says 12/12/2025, sample 10 says 12/12/2024.
# Pairing follows Sheet3 regardless of those discrepancies.
# ============================================================

# ----------------------------
# 0. PACKAGES
# ----------------------------
cran_packages <- c('dplyr', 'ggplot2', 'pheatmap', 'jsonlite', 'tibble')
bioc_packages <- c('tximport', 'DESeq2')
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos='https://cloud.r-project.org')
if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask=FALSE, update=FALSE)
suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(dplyr)
  library(ggplot2); library(pheatmap); library(jsonlite)
})

# ----------------------------
# 1. DIRECTORIES
# ----------------------------
salmon_dir <- '/data/Succai_Human_GDC_Indole/Human_GDC_Indole_Salmon'
tx2gene_file <- file.path(salmon_dir, 'tx2gene.csv')
output_dir <- file.path(salmon_dir, 'DESeq2_Paired_All20')
qc_dir <- file.path(output_dir, 'QC')
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

# ----------------------------
# 2. SAMPLES FROM WORKBOOK SHEET 5
# NOTE: exact quantification folder names MUST match these strings.
# For Salmon outputs named differently, change ONLY salmon_sample_id.
# ----------------------------
metadata <- data.frame(
  sample_number = seq_len(20),
  salmon_sample_id = c(
    '1_GDC_indlps_11_7_2024',
    '2_GDC_lps_11_7_2024',
    '3_169936_inlps_11_25_2024',
    '4_169936_lps_11_25_2024',
    '5_001716_inlps_11_25_2024',
    '6_001716_lps_11_25_2024',
    '7_169968_inlps_11_25_2024',
    '8_169968_lps_11_25_2024',
    '9_001528_in_plus_lps_12_12_2025',
    '10_001528_lps_12_12_2024',
    '11_12_13_2024_No11',
    '12_12_13_2024_No12',
    '13_12_16_2024_No9',
    '14_12_16_2024_No10',
    '15_1_9_2025_No3',
    '16_1_9_2025_No4',
    '17_1_9_2025_No8',
    '18_1_9_2025_No9',
    '19_1_9_2025_No13',
    '20_1_9_2025_No14'
  ),
  donor_id = factor(rep(sprintf('Donor_%02d', 1:10), each=2)),
  treatment = factor(rep(c('inlps','lps'), times=10), levels=c('lps','inlps')),
  stringsAsFactors=FALSE
)
metadata$sample_id <- metadata$salmon_sample_id
rownames(metadata) <- metadata$sample_id
sample_ids <- metadata$sample_id

# Verify full design and pairing before any import.
stopifnot(nrow(metadata)==20L, !anyDuplicated(sample_ids),
          nlevels(metadata$donor_id)==10L,
          identical(as.integer(table(metadata$treatment)), c(10L,10L)))
pair_table <- table(metadata$donor_id, metadata$treatment)
if (!all(pair_table[, 'lps']==1L & pair_table[, 'inlps']==1L)) {
  stop('Each donor must have exactly one LPS and one Indole+LPS sample.')
}
write.csv(metadata[, c('sample_number','sample_id','donor_id','treatment')],
          file.path(output_dir,'Sample_Metadata.csv'), row.names=FALSE)
write.csv(as.data.frame.matrix(pair_table), file.path(qc_dir,'Pair_Design_Check.csv'))
print(metadata[, c('sample_number','sample_id','donor_id','treatment')])

files <- file.path(salmon_dir, sample_ids, 'quant.sf')
names(files) <- sample_ids
missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop('Missing Salmon quant.sf files; confirm folders/quantification first:\n',
       paste(missing_files, collapse='\n'))
}
if (!file.exists(tx2gene_file)) stop('tx2gene mapping not found: ',tx2gene_file)
cat('\nVerified ',length(files),' Salmon samples from 10 treatment pairs.\n', sep='')

# ----------------------------
# 3. SALMON MAPPING QC
# ----------------------------
salmon_qc <- dplyr::bind_rows(lapply(sample_ids, function(s) {
  path <- file.path(salmon_dir,s,'aux_info','meta_info.json')
  if (!file.exists(path)) {
    warning('Salmon metadata missing: ',s)
    return(data.frame(sample_id=s, percent_mapped=NA_real_,
                      num_processed=NA_real_, num_mapped=NA_real_))
  }
  info <- jsonlite::fromJSON(path)
  val <- function(name) if (is.null(info[[name]])) NA_real_ else as.numeric(info[[name]])[1]
  data.frame(sample_id=s, percent_mapped=val('percent_mapped'),
             num_processed=val('num_processed'), num_mapped=val('num_mapped'))
})) %>% left_join(metadata[,c('sample_id','donor_id','treatment')],by='sample_id')
write.csv(salmon_qc,file.path(qc_dir,'Salmon_Mapping_QC.csv'),row.names=FALSE)
print(salmon_qc)
map_plot <- ggplot(salmon_qc,aes(x=reorder(sample_id,percent_mapped),
                                 y=percent_mapped,fill=treatment)) +
  geom_col(na.rm=TRUE)+coord_flip()+theme_classic()+
  labs(title='Salmon Mapping Rates (20 samples)',x='Sample',y='Mapping rate (%)')
ggsave(file.path(qc_dir,'Salmon_Mapping_Rates.pdf'),map_plot,width=11,height=9)

# ----------------------------
# 4. TX2GENE QC
# ----------------------------
tx2gene <- read.csv(tx2gene_file, stringsAsFactors=FALSE)
if (!all(c('TXNAME','GENEID') %in% names(tx2gene)))
  stop('tx2gene.csv requires TXNAME and GENEID columns')
tx2gene <- tx2gene %>% select(TXNAME,GENEID) %>%
  filter(!is.na(TXNAME),!is.na(GENEID),TXNAME!='',GENEID!='') %>% distinct()
if (anyDuplicated(tx2gene$TXNAME)) stop('Ambiguous transcript-to-gene mapping')
annotation_qc <- bind_rows(lapply(sample_ids,function(s) {
  quant_ids <- read.delim(files[[s]],stringsAsFactors=FALSE)$Name
  matched <- sum(quant_ids %in% tx2gene$TXNAME)
  data.frame(sample_id=s,total_transcripts=length(quant_ids),
             matched_transcripts=matched,match_percent=100*matched/length(quant_ids))
}))
write.csv(annotation_qc,file.path(qc_dir,'Transcript_Annotation_QC.csv'),row.names=FALSE)
print(annotation_qc)
if (any(is.na(annotation_qc$match_percent)) || any(annotation_qc$match_percent<80))
  stop('Low or undefined transcript annotation matching; check transcript versions')

# ----------------------------
# 5. TXIMPORT
# ----------------------------
cat('\nImporting Salmon quantifications...\n')
txi <- tximport(files=files,type='salmon',tx2gene=tx2gene,countsFromAbundance='no')
stopifnot(identical(colnames(txi$counts),rownames(metadata)))
write.csv(txi$counts,file.path(output_dir,'Gene_Estimated_Counts.csv'))
cat('Imported genes:',nrow(txi$counts),'samples:',ncol(txi$counts),'\n')

# ----------------------------
# 6. LIBRARY QC
# ----------------------------
raw_counts <- txi$counts
library_qc <- data.frame(sample_id=colnames(raw_counts),
                         estimated_gene_counts=colSums(raw_counts),
                         detected_genes=colSums(raw_counts>0)) %>%
  left_join(metadata[,c('sample_id','donor_id','treatment')],by='sample_id')
write.csv(library_qc,file.path(qc_dir,'Library_Size_Detected_Genes_QC.csv'),row.names=FALSE)
for (which in c('estimated_gene_counts','detected_genes')) {
  p <- ggplot(library_qc,aes(x=reorder(sample_id,.data[[which]]),
                             y=.data[[which]],fill=treatment))+
    geom_col()+coord_flip()+theme_classic()+labs(x='Sample',y=which,title=which)
  name <- if (which=='estimated_gene_counts') 'Library_Sizes.pdf' else 'Detected_Genes.pdf'
  ggsave(file.path(qc_dir,name),p,width=11,height=9)
}

# ----------------------------
# 7. DESEQ2 OBJECT + LOW COUNT FILTER
# Original criterion: >=10 counts in >=3 samples.
# ----------------------------
dds <- DESeqDataSetFromTximport(txi=txi,colData=metadata,
                               design=~ donor_id + treatment)
genes_before <- nrow(dds)
keep <- rowSums(counts(dds)>=10)>=3
dds <- dds[keep,]
if (nrow(dds)<2L) stop('Insufficient genes after filtering')
filter_qc <- data.frame(genes_before_filtering=genes_before,
                        genes_after_filtering=nrow(dds),
                        genes_removed=genes_before-nrow(dds))
write.csv(filter_qc,file.path(qc_dir,'Gene_Filtering_QC.csv'),row.names=FALSE)
print(filter_qc)

# ----------------------------
# 8. NORMALIZATION / BLIND VST FOR QC
# ----------------------------
dds <- estimateSizeFactors(dds)
vsd_qc <- varianceStabilizingTransformation(dds,blind=TRUE)
norm_log <- log2(counts(dds,normalized=TRUE)+1)
pdf(file.path(qc_dir,'Normalized_Count_Distributions.pdf'),width=14,height=8)
boxplot(as.data.frame(norm_log),las=2,col='lightblue',outline=FALSE,cex.axis=0.65,
        ylab='Log2(normalized count + 1)',main='Normalized Count Distributions')
dev.off()

# ----------------------------
# 9. CORRELATION + CLUSTERING QC
# ----------------------------
sample_cor <- cor(assay(vsd_qc),method='pearson')
write.csv(sample_cor,file.path(qc_dir,'Sample_Correlation_Matrix.csv'))
ann <- metadata[colnames(sample_cor),c('donor_id','treatment'),drop=FALSE]
pdf(file.path(qc_dir,'Sample_Correlation_Heatmap.pdf'),width=13,height=12)
pheatmap(sample_cor,annotation_col=ann,annotation_row=ann,
         clustering_distance_rows=as.dist(1-sample_cor),
         clustering_distance_cols=as.dist(1-sample_cor),
         main='Sample-to-Sample Pearson Correlation')
dev.off()
pdf(file.path(qc_dir,'Sample_Hierarchical_Clustering.pdf'),width=14,height=8)
plot(hclust(dist(t(assay(vsd_qc)))),main='Hierarchical Clustering',
     xlab='',sub='',cex=0.7)
dev.off()

# ----------------------------
# 10. PCA QC (SAMPLES + PAIR CONNECTIONS)
# ----------------------------
pca_data <- plotPCA(vsd_qc,intgroup=c('treatment','donor_id'),returnData=TRUE)
percentVar <- round(100*attr(pca_data,'percentVar'))
pca_data$sample_id <- as.character(pca_data$name)
write.csv(pca_data,file.path(qc_dir,'PCA_Coordinates.csv'),row.names=FALSE)
axis_labs <- labs(x=paste0('PC1: ',percentVar[1],'%'),
                  y=paste0('PC2: ',percentVar[2],'%'))
pca_plot <- ggplot(pca_data,aes(PC1,PC2,color=treatment))+
  geom_point(size=3.5)+geom_text(aes(label=sample_id),vjust=-0.8,size=2.2,
                                check_overlap=TRUE)+
  theme_classic()+axis_labs+labs(title='RNA-seq PCA (all samples)')
ggsave(file.path(qc_dir,'PCA_Treatment_Donor.pdf'),pca_plot,width=12,height=9)
paired_pca <- ggplot(pca_data,aes(PC1,PC2))+
  geom_line(aes(group=donor_id),color='grey45',linewidth=0.7,linetype='dashed')+
  geom_point(aes(color=treatment),size=3.5)+
  geom_text(aes(label=sample_id),vjust=-0.8,size=2.2,check_overlap=TRUE)+
  theme_classic()+axis_labs+
  labs(title='Paired PCA: Indole + LPS vs LPS',
       subtitle='Dashed lines connect sample pairs defined in workbook Sheet3')
ggsave(file.path(qc_dir,'PCA_Paired_Donors.pdf'),paired_pca,width=12,height=9)

# ----------------------------
# 11. PAIRED DIFFERENTIAL EXPRESSION
# ----------------------------
cat('\nRunning paired DESeq2 (10 donor blocks)...\n')
dds <- DESeq(dds)
res <- results(dds,contrast=c('treatment','inlps','lps'),alpha=0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df %>% arrange(padj)
write.csv(res_df,file.path(output_dir,'DEG_INLPS_vs_LPS_All.csv'),row.names=FALSE)
sig <- res_df %>% filter(!is.na(padj),padj<0.05)
up <- sig %>% filter(log2FoldChange>0)
down <- sig %>% filter(log2FoldChange<0)
write.csv(sig,file.path(output_dir,'DEG_INLPS_vs_LPS_FDR005.csv'),row.names=FALSE)
write.csv(up,file.path(output_dir,'Upregulated_INLPS.csv'),row.names=FALSE)
write.csv(down,file.path(output_dir,'Downregulated_INLPS.csv'),row.names=FALSE)
cat('Significant:',nrow(sig),'up:',nrow(up),'down:',nrow(down),'\n')

# ----------------------------
# 12. COOK'S DISTANCE AND RESULTS QC
# ----------------------------
cook_assay <- assays(dds)[['cooks']]
if (!is.null(cook_assay)) {
  cook_max <- apply(cook_assay,2,function(x) if (all(is.na(x))) NA_real_ else max(x,na.rm=TRUE))
  cook_qc <- data.frame(sample_id=names(cook_max),max_cooks_distance=as.numeric(cook_max)) %>%
    left_join(metadata[,c('sample_id','donor_id','treatment')],by='sample_id')
  write.csv(cook_qc,file.path(qc_dir,'Cooks_Distance_QC.csv'),row.names=FALSE)
  p <- ggplot(cook_qc,aes(x=reorder(sample_id,max_cooks_distance),
                          y=max_cooks_distance,fill=treatment))+
    geom_col()+coord_flip()+theme_classic()+
    labs(title="Maximum Cook's Distance per Sample",x='Sample',y="Max Cook's distance")
  ggsave(file.path(qc_dir,'Cooks_Distance_QC.pdf'),p,width=11,height=9)
}
result_qc <- data.frame(total_tested_genes=nrow(res_df),
                        genes_with_NA_pvalue=sum(is.na(res_df$pvalue)),
                        genes_with_NA_padj=sum(is.na(res_df$padj)),
                        significant_FDR005=nrow(sig))
write.csv(result_qc,file.path(qc_dir,'DESeq2_Results_QC.csv'),row.names=FALSE)

# ----------------------------
# 13. DISPERSION + MA PLOTS
# ----------------------------
pdf(file.path(qc_dir,'DESeq2_Dispersion_Plot.pdf'),width=9,height=7)
plotDispEsts(dds)
dev.off()
pdf(file.path(qc_dir,'DESeq2_MA_Plot.pdf'),width=9,height=7)
plotMA(res,ylim=c(-5,5),main='Indole + LPS vs LPS')
dev.off()

# ----------------------------
# 14. NORMALIZED COUNTS
# ----------------------------
write.csv(counts(dds,normalized=TRUE),file.path(output_dir,'Normalized_Gene_Counts.csv'))

# ----------------------------
# 15. VOLCANO
# ----------------------------
volcano_df <- res_df %>% filter(!is.na(padj),!is.na(log2FoldChange)) %>%
  mutate(significance=case_when(padj<0.05 & log2FoldChange>1 ~ 'Upregulated',
                                padj<0.05 & log2FoldChange< -1 ~ 'Downregulated',
                                TRUE ~ 'Not Significant'),
         neg_log10_padj=-log10(pmax(padj,.Machine$double.xmin)))
volcano <- ggplot(volcano_df,aes(log2FoldChange,neg_log10_padj,color=significance))+
  geom_point(alpha=0.65,size=1.4)+
  scale_color_manual(values=c(Upregulated='red',Downregulated='blue',
                              'Not Significant'='grey75'))+
  geom_hline(yintercept=-log10(0.05),linetype='dashed')+
  geom_vline(xintercept=c(-1,1),linetype='dashed')+theme_classic()+
  labs(title='Indole + LPS vs LPS (10 donor pairs)',
       subtitle='Highlight: FDR < 0.05 and |log2FC| > 1',
       x='Log2 Fold Change',y='-Log10 adjusted P-value')
ggsave(file.path(output_dir,'Volcano_INLPS_vs_LPS.pdf'),volcano,width=9,height=7)

# ----------------------------
# 16. TOP 50 DEG HEATMAP
# ----------------------------
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
top_genes <- head(sig$gene_id,50)
if (length(top_genes)>=2) {
  heatmap_mat <- assay(vsd)[top_genes,,drop=FALSE]
  ann <- metadata[colnames(heatmap_mat),c('donor_id','treatment'),drop=FALSE]
  pdf(file.path(output_dir,'Top50_DEG_Heatmap.pdf'),width=12,height=13)
  pheatmap(heatmap_mat,scale='row',annotation_col=ann,
           cluster_rows=TRUE,cluster_cols=TRUE,fontsize_row=7,
           main='Top 50 DEGs: Indole + LPS vs LPS')
  dev.off()
}

# ----------------------------
# 17. SAVE OBJECTS AND RUN DETAILS
# ----------------------------
saveRDS(dds,file.path(output_dir,'DESeq2_dds.rds'))
saveRDS(vsd,file.path(output_dir,'DESeq2_vst.rds'))
saveRDS(vsd_qc,file.path(qc_dir,'DESeq2_vst_blind_QC.rds'))
writeLines(capture.output(sessionInfo()),file.path(output_dir,'R_SessionInfo.txt'))
cat('\n===========================================\n',
    'ANALYSIS COMPLETE\n',
    'Samples: ',ncol(dds),'\n',
    'Donor pairs: ',nlevels(dds$donor_id),'\n',
    'Genes retained: ',nrow(dds),'\n',
    'Significant DEGs: ',nrow(sig),'\n',
    'Up: ',nrow(up),' | Down: ',nrow(down),'\n',
    'Results: ',output_dir,'\n',
    'QC: ',qc_dir,'\n',
    '===========================================\n',sep='')
