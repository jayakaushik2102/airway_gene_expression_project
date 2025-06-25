install.packages("tidyverse")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("EnhancedVolcano")
install.packages("RColorBrewer")
install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2",force = TRUE)
BiocManager::install("airway")
BiocManager::install("EnhancedVolcano")


library(tidyverse)
library(DESeq2)
library(airway)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

data("airway")
airway$dex <- relevel(airway$dex, "untrt")

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)
summary(res)

res_ordered <- res[order(res$pvalue), ]
head(res_ordered, 10)

dir.create("airway_gene_expression_project")

write.csv(as.data.frame(res), "airway_gene_expression_project/DESeq2_results.csv")

png("airway_gene_expression_project/volcano_plot.png", width = 1000, height = 800)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot - Airway DE Genes',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                col = c("grey", "blue", "orange", "red"))
dev.off()

vsd <- vst(dds, blind = FALSE)
topGenes <- head(order(res$pvalue), 20)
mat <- assay(vsd)[topGenes, ]

custom_colors <- colorRampPalette(c("white", "skyblue", "darkblue"))(100)

png("airway_gene_expression_project/heatmap_top_genes.png", width = 800, height = 600)
pheatmap(mat,
         scale = "row",
         show_rownames = TRUE,
         main = "Top 20 Differentially Expressed Genes",
         color = custom_colors)
dev.off()


