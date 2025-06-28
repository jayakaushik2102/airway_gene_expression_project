# Install essential packages for analysis and visualization
install.packages("tidyverse")        # For data manipulation and visualization
install.packages("pheatmap")         # For creating heatmaps
install.packages("ggplot2")          # For general plotting (part of tidyverse)
install.packages("EnhancedVolcano")  # For creating volcano plots
install.packages("RColorBrewer")     # For color palettes
install.packages("BiocManager")      # To install Bioconductor packages

# Check if BiocManager is installed; if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")     # Safe check for BiocManager install

# Install Bioconductor packages
BiocManager::install("DESeq2",force = TRUE)# For differential expression analysis
BiocManager::install("airway")# Sample dataset used in tutorials
BiocManager::install("EnhancedVolcano")# For better volcano plot visuals


library(tidyverse)# Load tidyverse for data wrangling
library(DESeq2)# Load DESeq2 for RNA-seq analysis
library(airway)# Load airway dataset
library(EnhancedVolcano)# Load package for volcano plot
library(pheatmap)# Load package for heatmap
library(RColorBrewer)# Load color palette support

# Load the sample data
data("airway")   # Load airway dataset
airway$dex <- relevel(airway$dex, "untrt")  # Relevel to set 'untrt' as the reference condition

# Create DESeq2 dataset and run differential expression analysis
dds <- DESeqDataSet(airway, design = ~ cell + dex)# Create DESeq2 object with design formula
dds <- DESeq(dds)# Run the DESeq pipeline
res <- results(dds)# Extract results table
summary(res)# View summary of DE results

# Order results by p-value
res_ordered <- res[order(res$pvalue), ]# Sort genes by p-value (ascending)
head(res_ordered, 10)# View top 10 differentially expressed genes

# Create a directory to store output files
dir.create("airway_gene_expression_project")# Make a new folder for saving outputs

# Save DESeq2 results to CSV
write.csv(as.data.frame(res), "airway_gene_expression_project/DESeq2_results.csv")# Save results as CSV

# Generate and save volcano plot as PNG
png("airway_gene_expression_project/volcano_plot.png", width = 1000, height = 800)# Open PNG device
EnhancedVolcano(res,# Create volcano plot
                lab = rownames(res),# Gene names as labels
                x = 'log2FoldChange',# X-axis: fold change
                y = 'pvalue',# Y-axis: p-value
                title = 'Volcano Plot - Airway DE Genes',# Title of plot
                pCutoff = 0.05,# p-value threshold
                FCcutoff = 1.5,# Fold change threshold
                col = c("grey", "blue", "orange", "red"))# Color scheme
dev.off()# Close PNG device and save the plot

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)# Transform data for visualization

# Extract top 20 DE genes based on p-value
topGenes <- head(order(res$pvalue), 20)# Get indices of top 20 genes
mat <- assay(vsd)[topGenes, ]# Extract expression values for top genes

# Define custom color palette for heatmap
custom_colors <- colorRampPalette(c("white", "skyblue", "darkblue"))(100)# Create gradient color palette

# Generate and save heatmap
png("airway_gene_expression_project/heatmap_top_genes.png", width = 800, height = 600)# Open PNG device
pheatmap(mat,  # Create heatmap
         scale = "row",  # Scale each gene row
         show_rownames = TRUE,  # Show gene names
         main = "Top 20 Differentially Expressed Genes",  # Title
         color = custom_colors) # Use custom colors 
dev.off()   # Close PNG device and save the plot


