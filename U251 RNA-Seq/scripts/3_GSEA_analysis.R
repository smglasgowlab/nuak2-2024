library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(fgsea)
library(KEGGREST)
library(tools)

# Load RDS files
res.df <- readRDS("intermediate_data/1_res_df.rds")

# Create ranked gene list
gene.list <- data.frame("Gene ID" = res.df$symbol,
                        "Rank" = res.df$stat)

gene.list <- gene.list %>% 
  arrange(desc(Rank)) %>% 
  deframe()

# Perform GSEA using GO biological processes
GSEA.GO.BP <- gseGO(geneList = gene.list,
                           ont = "BP",
                           OrgDb = "org.Hs.eg.db",
                           keyType = "SYMBOL",
                           minGSSize = 50)

# Import custom gene signature for MES subtype
MES.signature.2010 <- scan("VERHAAK_GLIOBLASTOMA_MESENCHYMAL.v2024.1.Hs.grp", 
                           what = "character", 
                           quiet = TRUE)
MES.signature.2010 <- MES.signature.2010[4:length(MES.signature.2010)]

PN.signature.2010 <- scan("VERHAAK_GLIOBLASTOMA_PRONEURAL.v2024.1.Hs.grp",
                          what = "character",
                          quiet = TRUE)
PN.signature.2010 <- PN.signature.2010[4:length(PN.signature.2010)]

# Import custom gene signature for GO:0001837
EMT.signature <- scan("GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.v2024.1.Hs.grp", 
                      what = "character", 
                      quiet = TRUE)
EMT.signature <- EMT.signature[4:length(EMT.signature)]

# Perform GSEA with custom annotations
TERM2GENE <- data.frame(TERM = c(rep("MES", length(MES.signature.2010)), 
                                 rep("PN", length(PN.signature.2010)),
                                 rep("EMT", length(EMT.signature))),
                        NAME = c(MES.signature.2010, PN.signature.2010, EMT.signature))
GSEA.custom <- GSEA(gene.list, TERM2GENE = TERM2GENE, pvalueCutoff = 1)


# Add gene symbols to GO BP and KEGG results
GSEA.GO.BP.df <- data.frame(setReadable(GSEA.GO.BP, 'org.Hs.eg.db', 'ENTREZID'))
GSEA.custom.df <- data.frame(GSEA.custom@result)

# Save RDS files
saveRDS(GSEA.GO.BP, file = "intermediate_data/3_GSEA_GO_BP.rds")
saveRDS(GSEA.custom, file = "intermediate_data/3_GSEA_custom.rds")

# Save CSV
write.csv(GSEA.GO.BP.df, "../excel_results/GSEA_GO_BP.csv", row.names = FALSE)
write.csv(GSEA.custom.df, "../excel_results/GSEA_custom.csv", row.names = FALSE)
