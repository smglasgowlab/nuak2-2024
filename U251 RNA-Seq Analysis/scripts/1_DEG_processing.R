library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(cowplot)

# Color palette for all figures
palette.colors = rev(c(hcl.colors(9, palette="Purple-Green")))
                   
# Read in raw counts
df <- read.csv("raw_counts.csv")

# Convert column 1 to row names
counts <- df[,-1]
colnames(counts) <- c("U251_WT_1", "U251_WT_2", "U251_WT_3",
                      "U251_N2CR_1", "U251_N2CR_2", "U251_N2CR_3") 
rownames(counts) <- df[,1]

# Remove all genes with a row sum of less than 10 reads
counts <- counts[which(rowSums(counts) > 10),]

# Create metadata dataframe
id <- c("U251_WT_1", "U251_WT_2", "U251_WT_3", "U251_N2CR_1", "U251_N2CR_2", "U251_N2CR_3")
condition <- c("WT", "WT", "WT", "N2CR", "N2CR", "N2CR")
col.data <- data.frame(id, condition)
col.data$condition <- factor(col.data$condition, levels=c("WT", "N2CR"))

# Create DESeq object with WT vs. N2KO as design
dds <- DESeqDataSetFromMatrix(countData=counts, 
                                  colData=col.data,
                                  design=~condition)
dds <- DESeq(dds)

# Extract DESeq results contrasting N2KO (numerator in log2FC) against WT (denominator)
res <- results(dds, contrast=c("condition", "N2CR", "WT"))

# Convert results to dataframe
res.df <- data.frame(res@listData)
rownames(res.df) <- rownames(counts)

# Map ENTREZ gene ID and gene symbols, add symbols to counts df 
res.df$symbol <- mapIds(org.Hs.eg.db, keys=rownames(res.df), 
                        keytype="ENTREZID", column="SYMBOL")
counts$symbol <- res.df$symbol 
res.df <- res.df[!is.na(res.df$symbol), ]

# See if U251 direction of expression correlated with TCGA-GBM
tcga.res.df <- read.csv('GlioVis Data/GlioVis_high_vs_low_NUAK2.csv')
res.df$correlation <- -res.df$log2FoldChange / 
  tcga.res.df$logFC[match(res.df$symbol, tcga.res.df$X)]
res.df <- res.df %>% 
  mutate(correlation=case_when(
    correlation > 0 ~ "Correlation",
    correlation < 0 ~ "No Correlation",
    TRUE ~ NA_character_
  ))

# Set parameters for filtering DEGs
DEGs <- res.df[res.df$padj < 0.05 & 
                 abs(res.df$log2FoldChange) > 2 & 
                 res.df$baseMean > 10, ]
DEGs <- DEGs[!is.na(DEGs$log2FoldChange), ]

# Sort DEGs into upregulated and downregulated
DEGs.upregulated <- DEGs[DEGs$log2FoldChange > 2, ]
DEGs.downregulated <- DEGs[DEGs$log2FoldChange < -2, ]

# Normalize transcript counts and calculate Z-score 
dds.norm <- counts(dds, normalized =  TRUE)
dds.z <- t(apply(dds.norm, 1, scale))
colnames(dds.z) <- c("U251_WT_1", "U251_WT_2", "U251_WT_3", 
                     "U251_N2CR_1", "U251_N2CR_2", "U251_N2CR_3")
rownames(dds.z) <- counts$symbol
dds.z <- dds.z[!is.na(rownames(dds.z)), ]

# Save RDS Files
saveRDS(palette.colors, file="intermediate_data/1_palette_colors.rds")
saveRDS(counts, file="intermediate_data/1_counts.rds")
saveRDS(dds, file="intermediate_data/1_dds.rds")
saveRDS(res.df, file="intermediate_data/1_res_df.rds")
saveRDS(DEGs, file="intermediate_data/1_DEGs.rds")
saveRDS(DEGs.upregulated, file="intermediate_data/1_DEGs_upregulated.rds")
saveRDS(DEGs.downregulated, file="intermediate_data/1_DEGs_downregulated.rds")
saveRDS(dds.z, file="intermediate_data/1_dds_z.rds")

# Save CSV
write.csv(res.df, "../excel_results/res_df.csv", row.names=FALSE)
write.csv(DEGs, "../excel_results/DEGs.csv", row.names=FALSE)