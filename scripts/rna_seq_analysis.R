## Include libraries

library(here)
library(DESeq2)
library(gplots)

##Read in the csv

infile <- read.csv(here('data', 'GSE223024_deseq2_normalized_genes.csv.gz')) 

##Remove duplicate gene names

infile <- infile[!duplicated(infile$Gene_name), ]

##Create sample names for the 3 AR replicates and 3 control replicates

sample_names <- c('LN_AR_1.', 'LN_AR_2.', 'LN_AR_3.', 'LN_NC_1.', 'LN_NC_2.', 'LN_NC_3.') 

##Create the counts data frame for DESeq2 and round the numbers to integers

counts <- round(infile[,sample_names])
row.names(counts) <- infile$Gene_name

##Create the samples data frame for DESeq2 and group by AR vs Control

samples <- data.frame(
'group' = c('ar', 'ar', 'ar', 'ctrl', 'ctrl', 'ctrl'),
row.names = names(counts))

##Load DESeq2 library and run analysis using created data frames

dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ group)
dds <- DESeq(dds)
ar_vs_nc <- results(dds, contrast = c('group', 'ar', 'ctrl'))

##Create data frame of significant genes

sig_ar_vs_nc <- ar_vs_nc[
	!is.na(ar_vs_nc$padj) &
	ar_vs_nc$padj < 0.05 &
	abs(ar_vs_nc$log2FoldChange) >= 1,
]

##Save the data frame for later and save number of rows

write.table(sig_ar_vs_nc, 'colin_byrne_ar_vs_nc_deseq.txt')
ar_sig_gene_num <- nrow(sig_ar_vs_nc)

##Read in HOXB13 csv

hox_infile <- read.csv(here('data', 'GSE275928_Cuff_Gene_Counts.csv.gz'))

##Remove duplicate gene ids

hox_infile <- hox_infile[!duplicated(hox_infile$Gene_ID), ]

##Create sample names for 4 HOXB13 replicates and 4 control replicates

hox_sample_names <- c('HOXB13_rep1', 'HOXB13_rep2', 'HOXB13_rep3', 'HOXB13_rep4', 'NT_rep1', 'NT_rep2', 'NT_rep3', 'NT_rep4')

##Create the counts data frame for DESeq2 and round the numbers to integers

hox_counts <- round(hox_infile[,hox_sample_names])
row.names(hox_counts) <- hox_infile$Gene_ID

#Create the samples data frame for DESeq2 and group by Hox vs Control

hox_samples  <- data.frame(
'group' = c('hox', 'hox', 'hox', 'hox', 'ctrl', 'ctrl', 'ctrl', 'ctrl'),
row.names = names(hox_counts))

##Run analysis using created data frames

hox_dds <- DESeqDataSetFromMatrix(hox_counts, hox_samples, design = ~ group)
hox_dds <- DESeq(hox_dds)
hox_vs_nt <- results(hox_dds, contrast = c('group', 'hox', 'ctrl'))

##Create data frame of significant genes

sig_hox_vs_nt <- hox_vs_nt[
	!is.na(hox_vs_nt$padj) &
	hox_vs_nt$padj < 0.05 &
	abs(hox_vs_nt$log2FoldChange) >= 1,
]

##Save the data frame for later and save number of rows

write.table(sig_hox_vs_nt, 'colin_byrne_hox_vs_nt_deseq.txt')
hox_sig_gene_num <- nrow(sig_hox_vs_nt)

##Cast the results to data frames to facilitate merging

ar_vs_nc <- data.frame(ar_vs_nc)
hox_vs_nt <- data.frame(hox_vs_nt)

##Merge on gene name (parameter 0)

ar_hox <- merge(ar_vs_nc, hox_vs_nt, by = 0, suffixes = c('.AR', '.HOX'))
row.names(ar_hox) <- ar_hox[,1]
ar_hox <- ar_hox[,-1]

##Replace NA values with 0 and 1 for log2FoldChange and padj respectively

ar_hox$log2FoldChange.HOX[ is.na(ar_hox$log2FoldChange.HOX) ] <- 0
ar_hox$log2FoldChange.AR[ is.na(ar_hox$log2FoldChange.AR) ] <- 0
ar_hox$padj.HOX[ is.na(ar_hox$padj.HOX) ] <- 1
ar_hox$padj.AR[ is.na(ar_hox$padj.AR) ] <- 1

##Save merged data frame for later and save number of rows

write.table(ar_hox, 'colin_byrne_ar_hox_merged.txt')
merged_genes_num <- nrow(ar_hox)

##Create new columns containing genes downregulated in AR and HOX

ar_hox$down_AR <- ar_hox$padj.AR < 0.05 & ar_hox$log2FoldChange.AR < 0
ar_hox$down_HOX <- ar_hox$padj.HOX < 0.05 & ar_hox$log2FoldChange.HOX < 0

##Find number of genes downregulated in AR, HOX, and both

genes_down_ar <- sum(ar_hox$down_AR & !ar_hox$down_HOX, na.rm = TRUE)
genes_down_hox <- sum(ar_hox$down_HOX & !ar_hox$down_AR, na.rm = TRUE)
genes_down_both <- sum(ar_hox$down_HOX & ar_hox$down_AR, na.rm = TRUE)

##Get table of signicant genes based on padj value for heatmap

sig_ar_hox <- ar_hox[
    (ar_hox$padj.AR < 0.05 | ar_hox$padj.HOX < 0.05), 
]
sig_ar_hox_genes <- nrow(sig_ar_hox)

##Create heatmap
png(here("results", "heatmap.png"))
heatmap.2(as.matrix(sig_ar_hox[, c("log2FoldChange.HOX", "log2FoldChange.AR")]),
                                col = bluered(50),
                                breaks = 0.1*(-25:25),
                                symbreaks = TRUE,
                                symkey = TRUE,
                                trace = 'none',
                                labRow = FALSE,
                                density.info = 'none',
                                key.xlab = 'log2 Fold Change over Ctrl',
                                labCol = c('HOXB13 KD', 'AR KD'),
                                main = 'Colin Byrne')
dev.off()
##Save results to summary.txt

sink(here("results", "rna_seq_summary.txt"))

cat("Total rows in AR knockdown table:", ar_sig_gene_num, "\n")
cat("Total rows in HOXB13 knockdown table:", hox_sig_gene_num, "\n")
cat("Total rows in merged table:", merged_genes_num, "\n")
cat("AR-only downregulated:", genes_down_ar, "\n")
cat("HOXB13-only downregulated:", genes_down_hox, "\n")
cat("Downregulated in both:", genes_down_both, "\n")
cat("Downregulated genes based on padj:", sig_ar_hox_genes)

sink()

cat("Analysis complete. Results saved in results/ folder.\n")