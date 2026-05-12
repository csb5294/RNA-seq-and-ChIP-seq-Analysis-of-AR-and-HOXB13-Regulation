cat("Starting full pipeline...\n")

# Install required packages if missing
if (!require("here")) install.packages("here")
if (!require("DESeq2")) {
    install.packages("BiocManager")
    BiocManager::install("DESeq2")
}
if (!require("gplots")) install.packages("gplots")

# Create output folder if it doesn't exist
if (!dir.exists("output")) {
    dir.create("output")
}

# Run the RNA-seq script
source(here::here("scripts", "rna_seqA_analysis.R"))

cat("Pipeline complete! Check the output/ folder.\n")