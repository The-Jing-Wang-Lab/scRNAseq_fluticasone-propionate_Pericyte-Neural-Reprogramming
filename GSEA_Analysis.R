##### Packages #####

library(Seurat)
library(hdf5r)
library(future)
library(devtools)
library(dplyr)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

##### Data Loading #####


# setting to correct working directory 
setwd("C:/Users/apsre/Documents/University/Fourth Year/Fall Sem/Jing_Wang_Lab/scRNAseq/jiwang_single_cell_2025-12-19/wang_251219_OCM")

# folder that contains one subfolder per sample (Set7_con, Set7_diff, Set7_midC, Set7_midF)
directory <- "per_sample_outs"

# read in filtered (high-quality) cell × gene count matrices from 10x outputs
con  <- Read10X_h5(file.path(directory, "Set7_con",  "count", "sample_filtered_feature_bc_matrix.h5"))
diff <- Read10X_h5(file.path(directory, "Set7_diff", "count", "sample_filtered_feature_bc_matrix.h5"))
midC <- Read10X_h5(file.path(directory, "Set7_midC", "count", "sample_filtered_feature_bc_matrix.h5"))
midF <- Read10X_h5(file.path(directory, "Set7_midF", "count", "sample_filtered_feature_bc_matrix.h5"))

# convert raw count matrices into Seurat objects
con <- CreateSeuratObject(con,
                          min.cells = 3, # removes genes expressed in <3 cells
                          min.features = 200, # removes cells with <200 genes (prob dead)
                          project = "Set7_con")
con <- RenameCells(con, add.cell.id = "con") # add con_ to each cell tag
head(colnames(con))

diff <- CreateSeuratObject(diff, 
                           min.cells = 3, # removes genes expressed in <3 cells
                           min.features = 200, # removes cells with >200 genes (prob dead)
                           project = "Set7_diff")
diff <- RenameCells(diff, add.cell.id = "diff") # add diff_ to each cell tag
head(colnames(diff))

midC <- CreateSeuratObject(midC, 
                           min.cells = 3, # removes genes expressed in <3 cells
                           min.features = 200, # removes cells with >200 genes (prob dead)
                           project = "Set7_midC")
midC <- RenameCells(midC, add.cell.id = "midC") # add midC_ to each cell tag
head(colnames(midC))

midF <- CreateSeuratObject(midF, 
                           min.cells = 3, # removes genes expressed in <3 cells
                           min.features = 200, # removes cells with >200 genes (prob dead)
                           project = "Set7_midF")
midF <- RenameCells(midF, add.cell.id = "midF") # add midC_ to each cell tag
head(colnames(midF))

##### Pseudobulk #####


con$condition <- "CTRL"
diff$condition <- "DIFF"
midC$condition <- "MIDC"
midF$condition <- "MIDF"

con$sample <- "con_1"
diff$sample <- "diff_1"
midC$sample <- "midC_1"
midF$sample <- "midF_1"

# aggregate all 4 samples
combined <- merge(con, y = list(diff, midC, midF))

combined$condition <- combined$orig.ident

DefaultAssay(combined) <- "RNA"

combined$cluster <- Idents(seurat_set7)

library(Matrix)

combined <- JoinLayers(combined)
counts <- GetAssayData(combined, layer = "counts")



# Aggregated expression profiles by condition, pseudo-pseudo-bulk

# prepping metadata
combined$condition <- combined$orig.ident
combined$condition <- combined$condition[colnames(counts)]
# sanity check
stopifnot(length(combined$condition) == ncol(counts))


# split cells by condition
cells_by_cond <- split(colnames(counts), combined$condition)

# pseudobulk aggregation (sum RNA)
pb_list <- lapply(cells_by_cond, function(cells) {
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

pb <- do.call(cbind, pb_list)

# check
dim(pb) # 24681 rows, 4 columns
colnames(pb) # condition names

# normalize

library(edgeR)
packageVersion("edgeR")

pb <- round(pb)  # ensure integer counts

y <- DGEList(counts = pb)
y <- calcNormFactors(y)

logCPM <- edgeR::cpm(y, log = TRUE)

##### GO Analysis, upregulated #####


# run all gsea function
run_gsea <- function(logCPM, group1, group2) {
  
  # logFC
  logFC <- logCPM[, group1] - logCPM[, group2]
  
  gene_list <- logFC
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # convert to ENTREZ
  gene_df <- bitr(names(gene_list),
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
  # merge while preserving order
  gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]  # keep first mapping only
  
  gene_list <- gene_list[gene_df$SYMBOL]
  names(gene_list) <- gene_df$ENTREZID
  
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # GSEA
  gsea_BP <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db,
                   ont = "BP", minGSSize = 10, maxGSSize = 500,
                   pvalueCutoff = 0.1, verbose = FALSE)
  
  gsea_MF <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db,
                   ont = "MF", minGSSize = 10, maxGSSize = 500,
                   pvalueCutoff = 0.1, verbose = FALSE)
  
  gsea_CC <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db,
                   ont = "CC", minGSSize = 10, maxGSSize = 500,
                   pvalueCutoff = 0.1, verbose = FALSE)
  
  return(list(BP = gsea_BP, MF = gsea_MF, CC = gsea_CC))
}

# get top terms function 
get_top_terms <- function(gsea_obj, n = 10, p_cutoff = 0.05) {
  
  df <- as.data.frame(gsea_obj)
  
  df_top <- df %>%
    filter(p.adjust < p_cutoff) %>% # p<0.05 significance cutoff
    mutate(Count = sapply(core_enrichment, function(x) {
      length(strsplit(x, "/")[[1]])
    })) %>%
    arrange(desc(Count)) %>% # rank by count
    slice_head(n = n) # will take the top n terms (10)
  
  return(df_top)
}


# making plot function
plot_gsea <- function(df_top, title, p_cutoff = 0.05) {
  
  # filter for significant terms
  df_top <- df_top %>% filter(p.adjust < p_cutoff)
  
  # plot
  ggplot(df_top, aes(x = Count,
                     y = reorder(Description, Count),
                     fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "red", high = "blue") + # adds colour gradient based on p-val
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # bigger, bold, centered title
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    labs(title = title, # labels
         x = "Count",
         y = "",
         fill = "p.adjust")
}


# run analyses
res_midC_ctrl <- run_gsea(logCPM, "Set7_midC", "Set7_con")
res_diff_ctrl <- run_gsea(logCPM, "Set7_diff", "Set7_con")
res_diff_midC <- run_gsea(logCPM, "Set7_diff", "Set7_midC")
res_midF_ctrl <- run_gsea(logCPM, "Set7_midF", "Set7_con")

res_midF_midC <- run_gsea(logCPM, "Set7_midF", "Set7_midC")

##### Generate plots ####
# run corresponding analysis before plotting

## MIDC vs CTRL ##

# BP
df_top <- get_top_terms(res_midC_ctrl$BP, n = 10)
plot_gsea(df_top, "MidC vs Ctrl - GO Biological Process")


# MF
df_top <- get_top_terms(res_midC_ctrl$MF, n = 10)
plot_gsea(df_top, "MidC vs Ctrl - GO Molecular Function")

# CC
df_top <- get_top_terms(res_midC_ctrl$CC, n = 10)
plot_gsea(df_top, "MidC vs Ctrl - GO Cellular Component")


## DIFF vs CTRL ##

# BP
df_top <- get_top_terms(res_diff_ctrl$BP, n = 10)
plot_gsea(df_top, "Diff vs Ctrl - GO Biological Process")

# MF
df_top <- get_top_terms(res_diff_ctrl$MF, n = 10)
plot_gsea(df_top, "Diff vs Ctrl - GO Molecular Function")

# CC
df_top <- get_top_terms(res_diff_ctrl$CC, n = 10)
plot_gsea(df_top, "Diff vs Ctrl - GO Cellular Component")



## DIFF vs MIDC ##

# BP
df_top <- get_top_terms(res_diff_midC$BP, n = 10)
plot_gsea(df_top, "Diff vs MidC - GO Biological Process")

# MF
df_top <- get_top_terms(res_diff_midC$MF, n = 10)
plot_gsea(df_top, "Diff vs MidC - GO Molecular Function")

# CC
df_top <- get_top_terms(res_diff_midC$CC, n = 10)
plot_gsea(df_top, "Diff vs MidC - GO Cellular Component")


## MIDF vs CTRL ##

# BP
df_top <- get_top_terms(res_midF_ctrl$BP, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Biological Process")

# MF
df_top <- get_top_terms(res_midF_ctrl$MF, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Molecular Function")

# CC
df_top <- get_top_terms(res_midF_ctrl$CC, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Cellular Component")


## MIDF vs MIDC ##
# res_midF_midC

# BP
df_top <- get_top_terms(res_midF_midC$BP, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Biological Process")

# MF
df_top <- get_top_terms(res_midF_midC$MF, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Molecular Function")

# CC
df_top <- get_top_terms(res_midF_midC$CC, n = 10)
plot_gsea(df_top, "MidF vs Ctrl - GO Cellular Component")

