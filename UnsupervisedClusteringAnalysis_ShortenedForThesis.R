##### Packages #####

library(Seurat)
library(hdf5r)
library(future)
library(devtools)
library(dplyr)
library(writexl)
library(scCustomize)
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


# Number of cells and genes
preQC <- merge(
  con,
  y = list(diff, midC, midF),
  add.cell.ids = c("CON", "DIFF", "MIDC", "MIDF"),
  project = "Set7"
)

table(preQC$orig.ident)
# number of cells 

tapply(preQC$nFeature_RNA, preQC$orig.ident, mean)
# average number of genes per sample



# Filter cells one condition at a time

## Filter con

con$percent.mito <- PercentageFeatureSet(con, pattern = "^MT-")
head(con$percent.mito)
# % of genes whose names match the pattern
# pattern = "^mt-": starts with mt-
# adds a column to tbx18 for % of mitochondrial genes
con$percent.ribo <- PercentageFeatureSet(con, pattern ="^Rp[sl]")
# % of genes whose names match the pattern
# pattern = "^Rp[sl]": starts with mt-
# adds a column to tbx18 for % of ribosomal genes

count.max <- round(mean(con$nCount_RNA) + 2 * sd(con$nCount_RNA), digits = -2)
# In normally distributed data, 95% of values lie within ±2 SD
# Count.max will calculate the upper threshold
# We want to exclude cells with abnormally high UMI (RNA) counts
# Could be: 
# doublets (two cells captured together)
# overly large or mutated
# transcriptionally abnormal
# Ng2_naive$nCount_RNA: gives the number of UMI (RNA) counts
# round([value], digits=-2): rounds the value to the nearest hundred 
count.min <- round(mean(con$nCount_RNA) - 2 * sd(con$nCount_RNA), digits = -2)
# Lower threshold 
# Cells with less RNA counts than this may be dying cells or empty droplets
feat.max <- round(mean(con$nFeature_RNA) + 2 * sd(con$nFeature_RNA), digits = -2)
# Same, but threshold for genes
# upper threshold
# Ng2_naive$nFeature_RNA: gives the number of types of RNAs (genes)
feat.min <- round(mean(con$nFeature_RNA) - 2 * sd(con$nFeature_RNA), digits = -2)
# lower threshold

# Set minimum parameters to 0 if negative value
if (count.min < 0){
  count.min <- 0
} else {
  count.min <- count.min
}

if (feat.min < 0){
  feat.min <- 0
} else {
  feat.min <- feat.min
}

# Filter out the cells
con <- subset(
  con, 
  subset = nFeature_RNA > feat.min &
    nFeature_RNA < feat.max & 
    nCount_RNA < count.max & 
    nCount_RNA > count.min & 
    percent.mito < 12)
# subsetting Ng2 dataset to exclude abnormal detections
# percent.mito < 12: high mitochondrial gene percentagel, dying/stressed/damaged cells

VlnPlot(con, c("nCount_RNA", "nFeature_RNA","percent.mito"), pt.size=0.1, ncol = 3)
# violin plots for nCount_RNA (total reads per cell), nFeature_RNA (genes detected per cell), & percent.mito (mitochondrial %)
# violin plot: box plot with a kernel density plot to show the distribution of numeric data, helps to check if outliers have been taken out
hist(con$percent.mito, breaks =100)
# histogram of # of cells vs mitochondrial gene percentage
# breaks =100, splits the histogram x-axis into 100 bins
FeatureScatter(con, feature1 = "nCount_RNA", feature2 = "percent.mito")
# scatter plot
# x-axis = nCount_RNA
# y-axis = percent.mito
# each point is a cell
FeatureScatter(con, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# scatter plot
# x-axis = nCount_RNA
# y-axis = nFeature_RNA
# each point is a cell
hist(con$nFeature_RNA, breaks =50)
hist(con$nCount_RNA, breaks =50)


## Filter diff

diff$percent.mito <- PercentageFeatureSet(diff, pattern = "^MT-")
# % of genes whose names match the pattern
# pattern = "^mt-": starts with mt-
# adds a column to tbx18 for % of mitochondrial genes
diff$percent.ribo <- PercentageFeatureSet(diff, pattern ="^Rp[sl]")
# % of genes whose names match the pattern
# pattern = "^Rp[sl]": starts with mt-
# adds a column to tbx18 for % of ribosomal genes

count.max <- round(mean(diff$nCount_RNA) + 2 * sd(diff$nCount_RNA), digits = -2)
# In normally distributed data, 95% of values lie within ±2 SD
# Count.max will calculate the upper threshold
# We want to exclude cells with abnormally high UMI (RNA) counts
# Could be: 
# doublets (two cells captured together)
# overly large or mutated
# transcriptionally abnormal
# Ng2_naive$nCount_RNA: gives the number of UMI (RNA) counts
# round([value], digits=-2): rounds the value to the nearest hundred 
count.min <- round(mean(diff$nCount_RNA) - 2 * sd(diff$nCount_RNA), digits = -2)
# Lower threshold 
# Cells with less RNA counts than this may be dying cells or empty droplets
feat.max <- round(mean(diff$nFeature_RNA) + 2 * sd(diff$nFeature_RNA), digits = -2)
# Same, but threshold for genes
# upper threshold
# Ng2_naive$nFeature_RNA: gives the number of types of RNAs (genes)
feat.min <- round(mean(diff$nFeature_RNA) - 2 * sd(diff$nFeature_RNA), digits = -2)
# lower threshold

# Set minimum parameters to 0 if negative value
if (count.min < 0){
  count.min <- 0
} else {
  count.min <- count.min
}

if (feat.min < 0){
  feat.min <- 0
} else {
  feat.min <- feat.min
}

# Filter out the cells
diff <- subset(
  diff, 
  subset = nFeature_RNA > feat.min &
    nFeature_RNA < feat.max & 
    nCount_RNA < count.max & 
    nCount_RNA > count.min & 
    percent.mito < 12)
# subsetting Ng2 dataset to exclude abnormal detections
# percent.mito < 12: high mitochondrial gene percentagel, dying/stressed/damaged cells

VlnPlot(diff, c("nCount_RNA", "nFeature_RNA","percent.mito"), pt.size=0.1, ncol = 3)
# violin plots for nCount_RNA (total reads per cell), nFeature_RNA (genes detected per cell), & percent.mito (mitochondrial %)
# violin plot: box plot with a kernel density plot to show the distribution of numeric data, helps to check if outliers have been taken out
hist(diff$percent.mito, breaks =100)
# histogram of # of cells vs mitochondrial gene percentage
# breaks =100, splits the histogram x-axis into 100 bins
FeatureScatter(diff, feature1 = "nCount_RNA", feature2 = "percent.mito")
# scatter plot
# x-axis = nCount_RNA
# y-axis = percent.mito
# each point is a cell
FeatureScatter(diff, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# scatter plot
# x-axis = nCount_RNA
# y-axis = nFeature_RNA
# each point is a cell
hist(diff$nFeature_RNA, breaks =50)
hist(diff$nCount_RNA, breaks =50)

## Filter midC

midC$percent.mito <- PercentageFeatureSet(midC, pattern = "^MT-")
# % of genes whose names match the pattern
# pattern = "^mt-": starts with mt-
# adds a column to tbx18 for % of mitochondrial genes
midC$percent.ribo <- PercentageFeatureSet(midC, pattern ="^Rp[sl]")
# % of genes whose names match the pattern
# pattern = "^Rp[sl]": starts with mt-
# adds a column to tbx18 for % of ribosomal genes

count.max <- round(mean(midC$nCount_RNA) + 2 * sd(midC$nCount_RNA), digits = -2)
# In normally distributed data, 95% of values lie within ±2 SD
# Count.max will calculate the upper threshold
# We want to exclude cells with abnormally high UMI (RNA) counts
# Could be: 
# doublets (two cells captured together)
# overly large or mutated
# transcriptionally abnormal
# Ng2_naive$nCount_RNA: gives the number of UMI (RNA) counts
# round([value], digits=-2): rounds the value to the nearest hundred 
count.min <- round(mean(midC$nCount_RNA) - 2 * sd(midC$nCount_RNA), digits = -2)
# Lower threshold 
# Cells with less RNA counts than this may be dying cells or empty droplets
feat.max <- round(mean(midC$nFeature_RNA) + 2 * sd(midC$nFeature_RNA), digits = -2)
# Same, but threshold for genes
# upper threshold
# Ng2_naive$nFeature_RNA: gives the number of types of RNAs (genes)
feat.min <- round(mean(midC$nFeature_RNA) - 2 * sd(midC$nFeature_RNA), digits = -2)
# lower threshold

# Set minimum parameters to 0 if negative value
if (count.min < 0){
  count.min <- 0
} else {
  count.min <- count.min
}

if (feat.min < 0){
  feat.min <- 0
} else {
  feat.min <- feat.min
}

# Filter out the cells
midC <- subset(
  midC, 
  subset = nFeature_RNA > feat.min &
    nFeature_RNA < feat.max & 
    nCount_RNA < count.max & 
    nCount_RNA > count.min & 
    percent.mito < 12)
# subsetting Ng2 dataset to exclude abnormal detections
# percent.mito < 12: high mitochondrial gene percentagel, dying/stressed/damaged cells

VlnPlot(midC, c("nCount_RNA", "nFeature_RNA","percent.mito"), pt.size=0.1, ncol = 3)
# violin plots for nCount_RNA (total reads per cell), nFeature_RNA (genes detected per cell), & percent.mito (mitochondrial %)
# violin plot: box plot with a kernel density plot to show the distribution of numeric data, helps to check if outliers have been taken out
hist(midC$percent.mito, breaks =100)
# histogram of # of cells vs mitochondrial gene percentage
# breaks =100, splits the histogram x-axis into 100 bins
FeatureScatter(midC, feature1 = "nCount_RNA", feature2 = "percent.mito")
# scatter plot
# x-axis = nCount_RNA
# y-axis = percent.mito
# each point is a cell
FeatureScatter(midC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# scatter plot
# x-axis = nCount_RNA
# y-axis = nFeature_RNA
# each point is a cell
hist(midC$nFeature_RNA, breaks =50)
hist(midC$nCount_RNA, breaks =50)


## Filter midF

midF$percent.mito <- PercentageFeatureSet(midF, pattern = "^MT-")
# % of genes whose names match the pattern
# pattern = "^mt-": starts with mt-
# adds a column to tbx18 for % of mitochondrial genes
midF$percent.ribo <- PercentageFeatureSet(midF, pattern ="^Rp[sl]")
# % of genes whose names match the pattern
# pattern = "^Rp[sl]": starts with mt-
# adds a column to tbx18 for % of ribosomal genes

count.max <- round(mean(midF$nCount_RNA) + 2 * sd(midF$nCount_RNA), digits = -2)
# In normally distributed data, 95% of values lie within ±2 SD
# Count.max will calculate the upper threshold
# We want to exclude cells with abnormally high UMI (RNA) counts
# Could be: 
# doublets (two cells captured together)
# overly large or mutated
# transcriptionally abnormal
# Ng2_naive$nCount_RNA: gives the number of UMI (RNA) counts
# round([value], digits=-2): rounds the value to the nearest hundred 
count.min <- round(mean(midF$nCount_RNA) - 2 * sd(midF$nCount_RNA), digits = -2)
# Lower threshold 
# Cells with less RNA counts than this may be dying cells or empty droplets
feat.max <- round(mean(midF$nFeature_RNA) + 2 * sd(midF$nFeature_RNA), digits = -2)
# Same, but threshold for genes
# upper threshold
# Ng2_naive$nFeature_RNA: gives the number of types of RNAs (genes)
feat.min <- round(mean(midF$nFeature_RNA) - 2 * sd(midF$nFeature_RNA), digits = -2)
# lower threshold

# Set minimum parameters to 0 if negative value
if (count.min < 0){
  count.min <- 0
} else {
  count.min <- count.min
}

if (feat.min < 0){
  feat.min <- 0
} else {
  feat.min <- feat.min
}

# Filter out the cells
midF <- subset(
  midF, 
  subset = nFeature_RNA > feat.min &
    nFeature_RNA < feat.max & 
    nCount_RNA < count.max & 
    nCount_RNA > count.min & 
    percent.mito < 12)
# subsetting Ng2 dataset to exclude abnormal detections
# percent.mito < 12: high mitochondrial gene percentagel, dying/stressed/damaged cells

VlnPlot(midF, c("nCount_RNA", "nFeature_RNA","percent.mito"), pt.size=0.1, ncol = 3)
# violin plots for nCount_RNA (total reads per cell), nFeature_RNA (genes detected per cell), & percent.mito (mitochondrial %)
# violin plot: box plot with a kernel density plot to show the distribution of numeric data, helps to check if outliers have been taken out
hist(midF$percent.mito, breaks =100)
# histogram of # of cells vs mitochondrial gene percentage
# breaks =100, splits the histogram x-axis into 100 bins
FeatureScatter(midF, feature1 = "nCount_RNA", feature2 = "percent.mito")
# scatter plot
# x-axis = nCount_RNA
# y-axis = percent.mito
# each point is a cell
FeatureScatter(midF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# scatter plot
# x-axis = nCount_RNA
# y-axis = nFeature_RNA
# each point is a cell
hist(midF$nFeature_RNA, breaks =50)
hist(midF$nCount_RNA, breaks =50)

##### Normalization, Integration, PCA, UMAPs #####


# Integrate data to get rid of batch effect

obj.list <- list(con, diff, midC, midF)
# Seurat’s integration functions expect a list of objects, not separate variables
# list() combines the 4 Seurat objects into a list

library(future)
options(future.globals.maxSize = 10 * 1024^3)
# the next line is copies objects and uses too much memory 

obj.list <- lapply(obj.list, function(x) {
  SCTransform(x, verbose = FALSE)})
# Run SCTransform on each dataset
# lapply(): applies a function to every item in a list
# function(x) is defined as SCTransform
obj.list <- lapply(obj.list, function(x) {
  DefaultAssay(x) <- "SCT"
  x })
# set default assay to SCT
# lapply(): applies a function to every item in a list
# function(x), SCT data is not the DefaultAssay 

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
# GENES
# creating a variable 'features' that will select and store genes that will be used to compare datasets
# object.list = obj.list: calls the variable that contains the datasets I want to integrate
# nfeatures = 3000: selects the top 3000 most variable genes shared across datasets (to avoid noise)
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features)
# needed for findIntegrationAnchors
# Recomputes residuals specifically for the 3000 integration genes, Stores them in a format integration expects, Makes sure all datasets are comparable
# SCTransform predicts the expected expression, it compares it to the actual observed expression
# Residual = observed value − predicted value
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
# CELLS
# creating a variable 'anchors' that will select and store cells that are similar between datasets
# Anchor: A pair of cells (from different datasets) that look biologically similar.
# object.list = obj.list: calls the variable that contains the datasets I want to integrate
# normalization.method = "SCT": reporting the normalization method we used earlier
# anchor.features = features: reporting the feature genes we're using
# NOTE: Canonical Correlation Analysis (CCA) was used by deafult in this command (as part of Seurat v4)
seurat_set7 <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# now we integrate the data!
# defined anchorset and how we normalized


# Reducing dimensionality

seurat_set7 <- RunPCA(seurat_set7, verbose = FALSE)
# machine learning technique used to reduce the dimensionality of a dataset

ElbowPlot(seurat_set7, ndims = 40)
# Determine percent of variation associated with each PC
pct <- seurat_set7[["pca"]]@stdev / sum(seurat_set7[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
pcs <- min(co1, co2)
pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

library(ggplot2)

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()



seurat_set7 <- RunUMAP(seurat_set7, dims=1:14)
seurat_set7 <- FindNeighbors(seurat_set7, dims=1:14)

seurat_set7 <- FindClusters(seurat_set7, resolution = c(0.6, 0.8, 1, 1.2))
p1 <- DimPlot(seurat_set7, group.by = "integrated_snn_res.0.6", label = T)
p2 <- DimPlot(seurat_set7, group.by = "integrated_snn_res.0.8", label = T)
p3 <- DimPlot(seurat_set7, group.by = "integrated_snn_res.1", label = T)
p4 <- DimPlot(seurat_set7, group.by = "integrated_snn_res.1.2", label = T)

p1 + p2 + p3 + p4
# number clusters: p1=5 ; p2=6 ; p3=10 ; p4=13

seurat_set7 <- FindClusters(seurat_set7, resolution = 0.6)
p6 <- DimPlot(seurat_set7, group.by = "integrated_snn_res.0.6", label = T)
p6 # 12 clusters

Idents(seurat_set7) <- "integrated_snn_res.0.6"
DefaultAssay(seurat_set7) <- "integrated"

# umap by condition
DimPlot(seurat_set7, group.by = "orig.ident")


# Number of cells and genes - Post Quality Check Table
table(seurat_set7$orig.ident)
# number of cells 
tapply(seurat_set7$nFeature_RNA, seurat_set7$orig.ident, mean)
# average number of genes per sample


# make excel file with marker genes per cluster
cluster_markers_set7 <- FindAllMarkers(seurat_set7,
                                       logfc.threshold = 0.25,
                                       only.pos = T)
library(dplyr)
top200set7<-cluster_markers_set7 %>%
  group_by(cluster)%>%
  top_n(200,avg_log2FC)
write.csv(top200set7, file="../cluster_200markers_set7(res0.7).csv", quote=F)


## Feature plot figure for thesis 
library(patchwork)

genes <- c("NOTCH3", "PDGFRB", "ABCC9", "ZIC1", "CD248",
           "VIM", "S100A6", "MN1", "TLE4", "CYP1B1-AS1", "HS3ST4")

plots <- lapply(genes, function(g) {
  FeaturePlot_scCustom(seurat_set7, features = g, na_color = "lightgrey") +
    
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      plot.title = element_text(size = 12, hjust = 0.5),
      
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.6, "cm")
    ) +
    
    guides(color = guide_colorbar(
      barheight = 5,
      barwidth = 0.8,
      ticks = FALSE
    ))
})

# Add 1 blank plot to make 12 total (4x3 grid)
plots[[12]] <- ggplot() + theme_void()

# Force 4 columns which will automatically give 3 rows
final_plot <- wrap_plots(plots, ncol = 3)
final_plot

ggsave("featureplots_grid.png", final_plot,
       width = 15, height = 20, dpi = 300)


## umap by cell identity

seurat_set7$celltype <- "Pericyte 1"
seurat_set7$celltype[Idents(seurat_set7) %in% c("0", "1", "9", "10")] <- "Pericyte 2"
seurat_set7$celltype[Idents(seurat_set7) %in% c("2", "3", "7", "8")] <- "Progenitor"
seurat_set7$celltype[Idents(seurat_set7) %in% c("6")] <- "Cortical Neuron"

celltype_colours <- c("Cortical Neuron" = "#53917e", "Pericyte 1" = "#e4572e", "Pericyte 2" = "#ffc857", "Progenitor" = "#4464ad")

# Output umap
p_celltype <- DimPlot(
  seurat_set7,
  group.by = "celltype",
  label = TRUE,
  cols = celltype_colours
) + ggtitle("Cell Identity")
p_celltype

## bar graph of % in each condition of each cell identity
# cell counts
cell_counts <- table(seurat_set7$orig.ident, seurat_set7$celltype)
cell_counts

# percent calulation
cell_percent <- prop.table(cell_counts, margin = 1) * 100

# make dataframe
cell_df <- as.data.frame(cell_percent)
colnames(cell_df) <- c("Condition", "CellType", "Percent")

# plot
library(ggplot2)

# add numbers on top of graphs
ggplot(cell_df, aes(x = CellType, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(Percent, 1)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4) +
  theme_classic() +
  ylab("% of cells") +
  xlab("") +
  ylim(0, max(cell_df$Percent) + 5)




