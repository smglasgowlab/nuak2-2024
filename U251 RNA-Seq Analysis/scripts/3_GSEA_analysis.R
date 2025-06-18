library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(fgsea)
library(KEGGREST)
library(tools)
library(biomaRt)

# Load RDS files
res.df <- readRDS("intermediate_data/1_res_df.rds")
proliferation.signature <- readRDS("gene_signatures/GOBP_CELL_POPULATION_PROLIFERATION.rds")
hippo.signature <- readRDS("gene_signatures/Hippo_signaling_pathway.rds")

# Create ranked gene list
gene.list <- data.frame("Gene ID"=res.df$symbol,
                        "Rank"=res.df$stat)

gene.list <- gene.list %>% 
  arrange(desc(Rank)) %>% 
  deframe()

# Perform GSEA using GO biological processes
GSEA.GO.BP <- gseGO(geneList=gene.list,
                           ont="BP",
                           OrgDb="org.Hs.eg.db",
                           keyType="SYMBOL",
                           minGSSize=50)

# Create gene list for KEGG
gene.list.KEGG <- data.frame("Gene ID"=rownames(res.df),
                        "Rank"=res.df$stat)

gene.list.KEGG <- gene.list.KEGG %>% 
  arrange(desc(Rank)) %>% 
  deframe()

# KEGG enrichment of GSEA 
GSEA.KEGG <- gseKEGG(gene.list.KEGG,
                     organism="hsa",
                     maxGSSize=200,
                     minGSSize=20)

# Perform GSEA with custom annotations
TERM2GENE <- data.frame(TERM = c(rep("Cell Population Proliferation", length(proliferation.signature)),
                                 rep("Hippo", length(hippo.signature))),
                        NAME = c(proliferation.signature, hippo.signature))
GSEA.custom <- GSEA(gene.list, TERM2GENE = TERM2GENE, pvalueCutoff=1)

# Add gene symbols to GO BP and KEGG results
GSEA.GO.BP.df <- data.frame(setReadable(GSEA.GO.BP, 'org.Hs.eg.db', 'ENTREZID'))
GSEA.KEGG.df <- data.frame(setReadable(GSEA.KEGG, 'org.Hs.eg.db', 'ENTREZID'))
GSEA.custom.df <- data.frame(GSEA.custom@result)

# Save RDS files
saveRDS(GSEA.KEGG, file="intermediate_data/3_GSEA_KEGG.rds")
saveRDS(GSEA.GO.BP, file="intermediate_data/3_GSEA_GO_BP.rds")
saveRDS(GSEA.custom, file="intermediate_data/3_GSEA_custom.rds")

# Save CSV
write.csv(GSEA.KEGG.df, "../excel_results/GSEA_KEGG.csv", row.names=FALSE)
write.csv(GSEA.GO.BP.df, "../excel_results/GSEA_GO_BP.csv", row.names=FALSE)
write.csv(GSEA.custom.df, "../excel_results/GSEA_custom.csv", row.names=FALSE)