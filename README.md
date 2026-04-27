# RNA-seq-and-ChIP-seq-Analysis-of-AR-and-HOXB13-Regulation

## Overview
This project analyzes RNA-seq and ChIP-seq data to investigate whether the transcription factor HOXB13 upstream of androgen receptor (AR) signaling in prostate cancer cells. Using public data sets from GEO, I compare expression and enhancer activity changes following AR and HOXB13 knockdown. All datasets were obtained from NCBI GEO. Raw or processed data were not redistributed in this repository due to size constraints; instead, direct download links are provided.

## Relation to Prior Work

This analysis is motivated by [Xiang et al 2025, CRISPR screening identifies regulators of enhancer-
mediated androgen receptor transcription in advanced prostate cancer, Cell Reports.](https://doi.org/10.1016/j.celrep.2025.115312) The paper proposes that HOXB13 functions as a master regulator upstream of AR in prostate cancer cells.

While the original study demonstrates that HOXB13 regulates AR expression and enhancer activity, it does not explicitly distinguish direct HOXB13 effects from indirect effects mediated through AR. This projects addresses that through analyzing differential expression between AR and HOXB13 knockdown, and analyzing overlap in enchancer activity changes using ChIP-seq data.

This approach allows for an independent evaluation of whether HOXB13-driven transcriptional changes are largely explained by downstream AR signaling.

## Biological Questions
If HOXB13 is a master regulator upstream of AR, then:
* Genes downregulated by AR knockdown should also be downregulated by HOXB13 knockdown
* AR-regulated enhancers should lose activity when HOXB13 is depleted

## Methods

### RNA-seq
* Differential expression analyzed in R using DESeq2
* Filtered using padj < 0.05 and log2FoldChange thresholds

### ChIP-Seq
* Peak intersection using BEDTools
* Defined enhancers by filtering out promoters
* Visualized using deepTools

## Results

## Data Sources

### RNA-seq
* AR knockdown (LNCaP cells)
  * GEO accession: [GSE223024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223024)
  * File used: GSE223024_deseq2_normalized_genes.csv.gz
* HOXB13 knockdown (LNCaP cells)
  * GEO accession: [GSE275928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275928)
  * File used: GSE275928_Cuff_Gene_Counts.csv.gz
### ChIP-seq
* ATAC-seq (open chromatin regions, LNCaP control)
  * GEO accession:[GSM7893414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7893414)
  * File used: GSM7893414_LNCaP_Vehicle-CFS-0_Rep3_TR_peaks.narrowPeak.gz
