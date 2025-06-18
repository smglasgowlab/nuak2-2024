library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(tools)

# Color palette for figures
palette.colors = c(hcl.colors(9, palette="Red-Blue"))

# Read in DEGs for NUAK2 high vs. low
res.df <- read.csv("GlioVis Data/GlioVis_high_vs_low_NUAK2.csv")

# Format res.df for upload
res.df.2 <- res.df[, c(3, 2, 4, 5, 6, 7, 1)]
colnames(res.df.2) <- c("baseMean", "log2FoldChange", "stat", 
                        "pvalue", "padj", "B", "symbol")

# Set first column to rownames
rownames(res.df) <- res.df$X
res.df <- res.df[,-1]

# Create an ENTREZID column
res.df$entrez <- mapIds(org.Hs.eg.db, 
                        keys = rownames(res.df), 
                        keytype="SYMBOL", 
                        column="ENTREZID",
                        multiVals="first")

# DEGs are LFC > 1 and p.adj < 0.05
DEGs <- res.df[abs(res.df$logFC) > 1 & res.df$adj.P.Val < 0.05, ]
DEGs$symbol <- rownames(DEGs)

# Separate into upregulated and downregualted DEGs
DEGs.upregulated <- DEGs[DEGs$logFC > 0, ]
DEGs.downregulated <- DEGs[DEGs$logFC < 0, ]

# Format DEGs for upload
DEGs.2 <- DEGs[, c(2, 1, 3, 4, 5, 6, 8)]
colnames(DEGs.2) <- c("baseMean", "log2FoldChange", "stat", 
                        "pvalue", "padj", "B", "symbol")

# Save RDS
saveRDS(palette.colors, file="intermediate_data/1_palette_colors.rds")
saveRDS(res.df, file="intermediate_data/1_res_df.rds")
saveRDS(DEGs, file="intermediate_data/1_DEGs.rds")
saveRDS(DEGs.upregulated, file="intermediate_data/1_DEGs_upregulated.rds")
saveRDS(DEGs.downregulated, file="intermediate_data/1_DEGs_downregulated.rds")

# Save CSV
write.csv(res.df.2, "../excel_results/res_df.csv", row.names = FALSE)
write.csv(DEGs.2, "../excel_results/DEGs.csv", row.names = FALSE)
