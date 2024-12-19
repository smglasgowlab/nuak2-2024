library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(tools)

# Load RDS files
counts <- readRDS("intermediate_data/1_counts.rds")
res.df <- readRDS("intermediate_data/1_res_df.rds")
DEGs <- readRDS("intermediate_data/1_DEGs.rds")
palette.colors <- readRDS("intermediate_data/1_palette_colors.rds")

# ORA for GO categories
GO.BP <- enrichGO(gene = rownames(DEGs),
                  universe = rownames(counts),
                  OrgDb = 'org.Hs.eg.db',
                  ont = "BP",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  minGSSize = 50)

GO.BP <- setReadable(GO.BP, "org.Hs.eg.db")

GO.BP.simplify <- clusterProfiler::simplify(GO.BP,
                                            cutoff = 0.7,
                                            by = "p.adjust",
                                            select_fun = min,
                                            measure = "Wang",
                                            semData = NULL)

GO.BP.simplify <- dropGO(GO.BP.simplify, term = c("GO:0032835", "GO:0072006", 
                                                  "GO:0001656", "GO:0007423"))

GO.CC <- enrichGO(gene = rownames(DEGs),
                  universe = rownames(counts),
                  OrgDb = 'org.Hs.eg.db',
                  ont = "CC",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  minGSSize = 50)

GO.CC <- setReadable(GO.CC, "org.Hs.eg.db")

GO.CC.simplify <- clusterProfiler::simplify(GO.CC,
                                            cutoff = 0.7,
                                            by = "p.adjust",
                                            select_fun = min,
                                            measure = "Wang",
                                            semData = NULL)

GO.MF <- enrichGO(gene = rownames(DEGs),
                  universe = rownames(counts),
                  OrgDb = 'org.Hs.eg.db',
                  ont = "MF",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  minGSSize = 50)


# Save RDS files
saveRDS(GO.BP.simplify, file = "intermediate_data/2_GO_BP_simplify.rds")
saveRDS(GO.CC.simplify, file = "intermediate_data/2_GO_CC_simplify.rds")
saveRDS(GO.MF.simplify, file = "intermediate_data/2_GO_MF_simplify.rds")

# Save CSV
write.csv(GO.BP@result, file = "../excel_results/GO_BP.csv", row.names = FALSE)
write.csv(GO.CC@result, file = "../excel_results/GO_CC.csv", row.names = FALSE)
write.csv(GO.BP.simplify@result, file = "../excel_results/GO_BP_simplify.csv", 
          row.names = FALSE)
write.csv(GO.CC.simplify@result, file = "../excel_results/GO_CC_simplify.csv",
          row.names = FALSE)