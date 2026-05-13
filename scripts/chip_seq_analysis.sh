#!/bin/bash

set -euo pipefail

mkdir -p results intermediate data

AR=data/ar_peaks.narrowPeak
HOX=data/hox_peaks.narrowPeak
K27=data/k27ac_peaks.broadPeak
ATAC=data/GSM7893414_LNCaP_Vehicle-CFS-0_Rep3_TR_peaks.narrowPeak.gz

GENES=data/hg38.gencode.bed12
CHRS=data/hg38_chr_sizes_full.txt

##Unzip downloaded file
if [ ! -f intermediate/atac.narrowPeak ]; then
    gunzip -c "$ATAC" > intermediate/atac.narrowPeak
fi

##Create promoter region file as 2000 bp up and downstream of TSS
bedtools flank -s -l 2000 -r 0 -i "$GENES" -g "$CHRS" | bedtools slop -s -l 0 -r 2000 -i - -g "$CHRS" > intermediate/promoters.bed

##Get AR peaks that overlap with open ATAC-seq regions, are marked by H3K27ac, and do not overlap with HOXB13 peaks
bedtools intersect -u -a "$AR" -b intermediate/atac.narrowPeak | bedtools intersect -u -a - -b "$K27" | \
bedtools intersect -v -a - -b "$HOX" > intermediate/ar_active_peaks.bed

##Separate AR peaks into promoter and enhancer regions
bedtools intersect -u -a intermediate/ar_active_peaks.bed -b intermediate/promoters.bed > intermediate/ar_atac_peaks_promoter.bed
bedtools intersect -v -a intermediate/ar_active_peaks.bed -b intermediate/promoters.bed > intermediate/ar_atac_peaks_enhancer.bed

##Get HOX peaks that overlap with open ATAC-seq regions, are marked by H3K27ac, and do not overlap with AR peaks
bedtools intersect -u -a "$HOX" -b intermediate/atac.narrowPeak | bedtools intersect -u -a - -b "$K27" | \
bedtools intersect -v -a - -b "$AR" > intermediate/hox_active_peaks.bed

##Separate AR peaks into promoter and enhancer regions
bedtools intersect -u -a intermediate/hox_active_peaks.bed -b intermediate/promoters.bed > intermediate/hox_atac_promoter.bed
bedtools intersect -v -a intermediate/hox_active_peaks.bed -b intermediate/promoters.bed > intermediate/hox_atac_enhancer.bed

##Get ATAC-seq open regions that are marked by H3K27ac, not in promoter regions, and not overlapped by either AR peaks or HOXB13 peaks
bedtools intersect -u -a intermediate/atac.narrowPeak -b "$K27" | bedtools intersect -v -a - -b intermediate/promoters.bed | \
bedtools intersect -v -a - -b "$AR" | bedtools intersect -v -a - -b "$HOX" > intermediate/control_enhancer_set.bed

##Get the first 2000 entries as a subset
head -n 2000 intermediate/control_enhancer_set.bed > intermediate/control_enhancer_subset.bed

##Compute matrix
computeMatrix reference-point --referencePoint center \
-b 2000 -a 2000 \
--sortRegions descend \
--missingDataAsZero \
--binSize 50 \
-o results/k27.matrix \
-S data/ar_k27ac.bw data/hox_k27ac.bw data/ctrl_k27ac.bw \
-R intermediate/ar_atac_peaks_enhancer.bed intermediate/hox_atac_enhancer.bed intermediate/control_enhancer_set.bed

##Create heatmap
plotHeatmap -m results/k27.matrix -o results/k27_heatmap.png \
--zMax 1 \
--heatmapHeight 20 \
--colorMap Blues \
--whatToShow "heatmap and colorbar" \
--xAxisLabel "Peak dist (bp)" \
--samplesLabel "AR KD" "HOXB13 KD" Ctrl \
--regionsLabel "AR enh" "HOXB13 enh" "Ctrl enh" \
--plotTitle "Colin Byrne"

##Create profile plot
plotProfile -m results/k27.matrix -o results/k27_profile.png \
--plotWidth 8 \
--yAxisLabel "Avg H3K27ac CPM" \
--samplesLabel "AR KD" "HOXB13 KD" Ctrl \
--regionsLabel "AR enh" "HOXB13 enh" "Ctrl enh" \
--plotTitle "Colin Byrne"