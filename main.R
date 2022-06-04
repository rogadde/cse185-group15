library(Seurat)
library(tidyverse)
library(scCustomize)
library(SeuratDisk)

# Load in the matrices and combine them into a single object
data <- Read10X_GEO(data_dir = "C:/Users/rohin/Desktop/CSE 185/data")
merged <- Merge_Sparse_Data_All(matrix_list = data, add_cell_ids = names(data))
sample <- CreateSeuratObject(counts = merged, min.cells = 3, min.features = 200)

rm(data, merged) # remove objects from memory

# Add a column for % mitochondrial genes
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

# Filter the data based on the paper's parameters
sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize and variance-stabilize data
sample <- SCTransform(sample, conserve.memory = TRUE, vars.to.regress = "percent.mt")

# Annotate each condition
condition <- function(x) {
  if (x %in% c("GSM5746268", "GSM5746269")) {"Control"}
  else if (x %in% c("GSM5746263", "GSM5746265")) {"Dead_D0"}
  else if (x %in% c("GSM5746264", "GSM5746266")) {"Dead_D7"}
  else if (x %in% c("GSM5746259", "GSM5746261")) {"Alive_D0"}
  else if (x %in% c("GSM5746260", "GSM5746262")) {"Alive_D7"}
}

sample@meta.data$condition <- mapply(condition, sample@meta.data$orig.ident)

disease <- function(x) {
  if (x == "Control) {"Control"}
  else {"COVID-19"}
}

sample@meta.data$diseaseStatus <- mapply(disease, sample@meta.data$condition)

# Load in the reference and map the query onto the reference in order to get annotations
reference <- LoadH5Seurat("Downloads/pbmc_multimodal.h5seurat")

anchors <- FindTransferAnchors(
  reference = reference, 
  query = sample, 
  normalization.method = "SCT", 
  reference.reduction = "spca", 
  dims = 1:50
)

sample <- MapQuery(
  anchorset = anchors,
  query = sample,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

rm(reference) # remove reference from memory

# Add UMAP computations to sample
sample <- RunUMAP(sample, reduction = 'ref.spca', dims = 1:50)

# Postprocessing: subset sample such that there are no erythrocytes
Idents(object=sample) <- "predicted.celltype.l2"
sample <- subset(sample, idents = c(
  "CD14 Mono", "CD16 Mono", "CD4 TCM", "CD8 TEM", "MAIT", "CD4 CTL", "NK", "NK Proliferating",
  "Treg", "CD4 Naive", "B naive", "CD4 TEM", "gdT", "B intermediate", "pDC", "B memory", "NK_CD56bright", "CD8 Naive", "HSPC",
  "Plasmablast", "CD8 TCM", "CD4 Proliferating", "cDC1", "dnT", "cDC2", "ILC", "CD8 Proliferating", "ASDC", "Platelet")
)

# Generate UMAP plots for Fig. 1c and 1d
DimPlot(sample, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + labs(title = "")
DimPlot(sample, group.by = "predicted.celltype.l2", split.by = "diseaseStatus", label = FALSE, label.size = 3, repel = TRUE, ncol=1) + labs(title = "")

sample$condition <- factor(sample$condition, levels = c("Control", "Alive_D0", "Alive_D7", "Dead_D0", "Dead_D7"), ordered = TRUE)
DimPlot(sample, group.by = "predicted.celltype.l2", split.by = "condition", label = FALSE, label.size = 3, repel = TRUE) + labs(title = "")

# Conduct differential gene expression analysis
# Our analysis gets both upregulated and downregulated DEGs

# Fig. 2b - Control vs. Critical COVID-19 Day 7
disease <- function(x) {
  if (x == "Control) {"Control"}
  else if (x == "Alive_D0" | x == "Dead_D0") {"Critical_D0"}
  else if (x == "Alive_D7" | x == "Dead_D7") {"Critical_D7"}
}

sample@meta.data$diseaseStatus <- mapply(disease, sample@meta.data$condition)

# Split sample by celltype
Idents(sample) <- "diseaseStatus"
list <- SplitObject(sample, split.by="predicted.celltype.l1")

# Find DEGs
control_critical7 <- lapply(X = list, FUN = FindMarkers, ident.1 = "Control", ident.2 = "Critical_D7", test.use = "wilcox")

celltypes <- c("Mono", "B", "NK", "CD4 T", "CD8 T", "other", "DC", "other T")
range <- 1:length(celltypes)

# Label DEGs
for (i in range) { 
  name <- celltypes[i] 
  control_critical7[i][[name]]$cell <- name
  control_critical7[i][[name]]$sig <- ifelse(control_critical7[i][[name]]$p_val_adj < 0.05 & (control_critical7[i][[name]]$avg_log2FC > 0.58 | control_critical7[i][[name]]$avg_log2FC < -0.58), "Significant","Not Significant")
}

# Combine each list element into a data frame
data <- bind_rows(control_critical7)
df_count <- data %>% group_by(sig, cell) %>% count()
df_count <- data.frame(df_count)
df_count <- filter(df_count, df_count$sig == "Significant")

# Plot
ggplot(df_count, aes(x = reorder(cell, -n), y = n)) + 
  geom_col(fill="red3") + 
  labs(x = "Cell Type", y = "Number of Genes", title = "Control vs. Critical COVID-19 Day 7") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Fig. 2d
Idents(sample) <- "condition"
list <- SplitObject(sample, split.by="predicted.celltype.l1")
a7_d7 <- lapply(X = list, FUN = FindMarkers, ident.1 = "Alive_D7", ident.2 = "Dead_D7", test.use = "wilcox")

for (i in range) { 
  name <- celltypes[i] 
  a7_d7[i][[name]]$cell <- name
  a7_d7[i][[name]]$sig <- ifelse(a7_d7[i][[name]]$p_val_adj < 0.05 & (a7_d7[i][[name]]$avg_log2FC > 0.58 | a7_d7[i][[name]]$avg_log2FC < -0.58), "Significant","Not Significant")
}

data <- bind_rows(a7_d7)
df_count <- data %>% group_by(sig, cell) %>% count()
df_count <- data.frame(df_count)
df_count <- filter(df_count, df_count$sig == "Significant")

ggplot(df_count, aes(x = reorder(cell, -n), y = n)) + 
  geom_col(fill="red3") + 
  labs(x = "Cell Type", y = "Number of Genes", title = "Alive vs. Deceased COVID-19 Day 7") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
# Fig. 3a
# Get B cell subsets for 3 conditions
Idents(sample) <- "predicted.celltype.l2"
currCell <- subset(sample, idents = c("B naive", "B intermediate", "B memory", "Plasmablast"))
Idents(currCell) <- "condition"
currCell_subset <- subset(currCell, idents = c("Control", "Alive_D7", "Dead_D7"))
currCell_subset$condition <- factor(currCell_subset$condition, levels = c("Control", "Alive_D7", "Dead_D7"), ordered = TRUE)

# Subset by condition
Idents(currCell_subset) <- "condition"
control <- subset(currCell_subset, idents = c("Control"))
alive7 <- subset(currCell_subset, idents = c("Alive_D7"))
dead7 <- subset(currCell_subset, idents = c("Dead_D7"))

# Generate UMAPs
DimPlot(control, group.by = "predicted.celltype.l2") + labs(title = "Control") + NoLegend()
DimPlot(alive7, group.by = "predicted.celltype.l2") + labs(title = "Alive Day 7") + NoLegend()
DimPlot(dead7, group.by = "predicted.celltype.l2") + xlim(-12, 2) + labs(title = "Dead Day 7")
