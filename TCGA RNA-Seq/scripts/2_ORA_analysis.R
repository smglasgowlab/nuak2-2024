library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(tools)

# Load RDS files
res.df <- readRDS("intermediate_data/1_res_df.rds")
DEGs <- readRDS("intermediate_data/1_DEGs.rds")

# ORA for GO:BP and GO:CC
GO.BP <- enrichGO(gene = rownames(DEGs),
                  universe = rownames(res.df),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = "SYMBOL",
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

GO.CC <- enrichGO(gene = rownames(DEGs),
                  universe = rownames(res.df),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = "SYMBOL",
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

# Save RDS 
saveRDS(GO.BP.simplify, file = "intermediate_data/2_GO_BP_simplify.rds")
saveRDS(GO.CC.simplify, file = "intermediate_data/2_GO_CC_simplify.rds")

# Save CSV
write.csv(GO.BP@result, file = "../excel_results/GO_BP.csv", row.names = FALSE)
write.csv(GO.CC@result, file = "../excel_results/GO_CC.csv", row.names = FALSE)
write.csv(GO.BP.simplify@result, file = "../excel_results/GO_BP_simplify.csv", 
          row.names = FALSE)
write.csv(GO.CC.simplify@result, file = "../excel_results/GO_CC_simplify.csv",
          row.names = FALSE)