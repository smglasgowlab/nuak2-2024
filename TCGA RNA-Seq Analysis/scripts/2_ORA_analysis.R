library(tidyverse)
library(conflicted)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(tools)

# Load RDS files
res.df <- readRDS("intermediate_data/1_res_df.rds")
DEGs.upregulated <- readRDS("intermediate_data/1_DEGs_upregulated.rds")
DEGs.downregulated <- readRDS("intermediate_data/1_DEGs_downregulated.rds")

GO_analysis <- function(ont, direction) {
  GO <- paste0("GO.", ont, ".", direction)
  DEGs.df <- paste0("DEGs.", direction)
  
  # GO analysis 
  assign(GO, 
    enrichGO(gene=rownames(get(DEGs.df)), 
             universe=rownames(res.df), 
             OrgDb='org.Hs.eg.db', 
             keyType="SYMBOL", 
             ont=ont, 
             pvalueCutoff =0.01, 
             qvalueCutoff=0.05, 
             minGSSize=50), 
    envir = .GlobalEnv
  )
  
  # Simplify GO terms similarity with score 0.7
  assign(paste0(GO, ".simplify"), 
         clusterProfiler::simplify(get(GO),
                                   cutoff=0.7,
                                   by="p.adjust",
                                   select_fun=min,
                                   measure="Wang",
                                   semData=NULL),
         envir = .GlobalEnv
  )
  
  # Save RDS, csv files 
  GO.name <- paste0("GO_", ont, "_", direction)
  GO.simplify.name <- paste0("GO_", ont, "_", direction, "_simplify")
  
  saveRDS(get(GO),
          file=paste0("intermediate_data/2_GO_", ont, "_", direction, ".rds"))
  saveRDS(get(paste0(GO, ".simplify")),
          file=paste0("intermediate_data/2_GO_", ont, "_", direction, "_simplify.rds"))
  
  write.csv(get(GO)@result,
            file=paste0("../excel_results/", GO.name, ".csv"))
  write.csv(get(paste0(GO, ".simplify"))@result, 
            file=paste0("../excel_results/", GO.simplify.name, ".csv"))
}

# GO analysis for BP, MF ontologies, set downregulated ratios to negative and merge
directions = c("upregulated", "downregulated")
onts = c("BP", "CC")
for (ont in onts) {
  for (direction in directions) {
    GO_analysis(ont, direction)
  }
}