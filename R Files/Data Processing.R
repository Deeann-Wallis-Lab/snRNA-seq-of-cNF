data <- Read10X_h5(filename =  "/filtered_feature_bc_matrix.h5")

g1 <- CreateSeuratObject(counts = data$`Gene Expression`)
g1[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`)

DefaultAssay(g1) <- "RNA"

g1 <- NormalizeData(g1)
g1 <- FindVariableFeatures(g1, selection.method = "mean.var.plot")
g1 <- ScaleData(g1, features = VariableFeatures(g1))

DefaultAssay(g1) <- "HTO"
g1 <- NormalizeData(g1, assay = "HTO", normalization.method = "CLR")
g1 <- HTODemux(g1, assay = "HTO", positive.quantile = 0.99)

table(g1$HTO_classification.global)

Idents(g1) <- "HTO_maxID"
RidgePlot(g1, assay = "HTO", features = rownames(g1[["HTO"]])[1:3], ncol = 3)
Idents(g1) <- "HTO_classification.global"

g1.singlet <-  subset(g1, idents = "Singlet")

Idents(g1.singlet) <- "HTO_classification"

cNF15v1 <- subset(g1.singlet, idents = "15v1")
cNF15v1$orig.ident <- "cNF15v1"

cNF15v4 <- subset(g1.singlet, idents = "15v4")
cNF15v4$orig.ident <- "cNF15v4"


# next pooled sample

data2 <- Read10X_h5(filename =  "/filtered_feature_bc_matrix.h5")

g2 <- CreateSeuratObject(counts = data2$`Gene Expression`)
g2[["HTO"]] <- CreateAssayObject(counts = data2$`Antibody Capture`)

DefaultAssay(g2) <- "RNA"

g2 <- NormalizeData(g2)
g2 <- FindVariableFeatures(g2, selection.method = "mean.var.plot")
g2 <- ScaleData(g2, features = VariableFeatures(g2))

DefaultAssay(g2) <- "HTO"
g2 <- NormalizeData(g2, assay = "HTO", normalization.method = "CLR")
g2 <- HTODemux(g2, assay = "HTO", positive.quantile = 0.99)

table(g2$HTO_classification.global)

Idents(g2) <- "HTO_maxID"
RidgePlot(g2, assay = "HTO", features = rownames(g2[["HTO"]])[1:4], ncol = 4)
Idents(g2) <- "HTO_classification.global"

g2.singlet <-  subset(g2, idents = "Singlet")

Idents(g2.singlet) <- "HTO_classification"

cNF1v1 <- subset(g2.singlet, idents = "1v1")
cNF1v1$orig.ident <- "cNF1v1"

cNF1vc3 <- subset(g2.singlet, idents = "1vc3")
cNF1vc3$orig.ident <- "cNF1vc3"


# next pooled sample

data3 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g3 <- CreateSeuratObject(counts = data3$`Gene Expression`)
g3[["HTO"]] <- CreateAssayObject(counts = data3$`Antibody Capture`)

DefaultAssay(g3) <- "RNA"

g3 <- NormalizeData(g3)
g3 <- FindVariableFeatures(g3, selection.method = "mean.var.plot")
g3 <- ScaleData(g3, features = VariableFeatures(g3))

DefaultAssay(g3) <- "HTO"
g3 <- NormalizeData(g3, assay = "HTO", normalization.method = "CLR")
g3 <- HTODemux(g3, assay = "HTO", positive.quantile = 0.99)

table(g3$HTO_classification.global)

Idents(g3) <- "HTO_maxID"
RidgePlot(g3, assay = "HTO", features = rownames(g3[["HTO"]])[1:4], ncol = 4)
Idents(g3) <- "HTO_classification.global"

g3.singlet <-  subset(g3, idents = "Singlet")

Idents(g3.singlet) <- "HTO_classification"

cNF13v1 <- subset(g3.singlet, idents = "13v1")
cNF13v1$orig.ident <- "cNF13v1"

cNF13v4 <- subset(g3.singlet, idents = "13v4")
cNF13v4$orig.ident <- "cNF13v4"


# next pooled sample

data5 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g5 <- CreateSeuratObject(counts = data5$`Gene Expression`)
g5[["HTO"]] <- CreateAssayObject(counts = data5$`Antibody Capture`)

DefaultAssay(g5) <- "RNA"

g5 <- NormalizeData(g5)
g5 <- FindVariableFeatures(g5, selection.method = "mean.var.plot")
g5 <- ScaleData(g5, features = VariableFeatures(g5))

DefaultAssay(g5) <- "HTO"
g5 <- NormalizeData(g5, assay = "HTO", normalization.method = "CLR")
g5 <- HTODemux(g5, assay = "HTO", positive.quantile = 0.99)

table(g5$HTO_classification.global)

Idents(g5) <- "HTO_maxID"
RidgePlot(g5, assay = "HTO", features = rownames(g5[["HTO"]])[1:4], ncol = 4)
Idents(g5) <- "HTO_classification.global"

g5.singlet <-  subset(g5, idents = "Singlet")

Idents(g5.singlet) <- "HTO_classification"

cNF5v1 <- subset(g5.singlet, idents = "5v1")
cNF5v1$orig.ident <- "cNF5v1"

cNF5v4 <- subset(g5.singlet, idents = "5v4")
cNF5v4$orig.ident <- "cNF5v4"


# put groups together

treated <- merge(x = cNF13v4, y = c(cNF15v4, cNF1vc3, cNF5v4), add.cell.ids = c("cNF13v4", "cNF15v4", "cNF1vc3",  "cNF5v4"), project = "treated")


treated$orig.ident <- "treated"

control <- merge(x = cNF13v1, y = c(cNF15v1, cNF1v1,  cNF5v1 ), add.cell.ids = c("cNF13v1", "cNF15v1", "cNF1v1",  "cNF5v1"), project = "control")
control$orig.ident <- "control"


all_samples <- merge(x = treated, y = c(control), add.cell.ids = c("treated", "control"), project = "all")

#saveRDS(all_samples, file = "/samples_15_13_1_5_merged.RDS")


DefaultAssay(all_samples) <- "RNA"

# Quality control

all_samples[["percent.mt"]] <- PercentageFeatureSet(all_samples, pattern = "^MT-")

# Visualize QC metrics as a violin plot
DimPlot(harmony_all_samples_integrated_cell_type, split.by = "orig.ident")

#all_samples <- harmony_all_samples_integrated_cell_type

g1 <- VlnPlot(all_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 5, pt.size = 0)

ggsave(plot = g1,
       filename = "QC1.png",
       path = ".../samples_15_13_1_5/",
       width = 20,
       height = 5,
       dpi = 300, limitsize = F)
g1 <- VlnPlot(all_samples, features = "percent.mt", ncol = 5, pt.size = 0)
ggsave(plot = g1,
       filename = "percent.mt.png",
       path = ".../samples_15_13_1_5/",
       width = 20,
       height = 5,
       dpi = 300, limitsize = F)

g1 <- VlnPlot(all_samples, features = "nFeature_RNA", ncol = 5, pt.size = 0)
ggsave(plot = g1,
       filename = "nFeature_RNA.png",
       path = ".../samples_15_13_1_5/",
       width = 20,
       height = 5,
       dpi = 300, limitsize = F)

plot1 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

all_samples <- subset(all_samples, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# continue onto intagrating the data

Test.list <- SplitObject(all_samples, split.by = "orig.ident")

all_samples.list <- Test.list[c("treated", "control")]

# Normalize the data for Intergration step

for (i in 1:length(all_samples.list)) {
  all_samples.list[[i]] <- NormalizeData(all_samples.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  all_samples.list[[i]] <- FindVariableFeatures(all_samples.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}
# Create a list of Anchor points for intergration

reference.list <- all_samples.list[c("treated", "control")]

all_samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)

# intergrate the data
all_samples.integrated <- IntegrateData(anchorset = all_samples.anchors, dims = 1:50)

DefaultAssay(all_samples.integrated) <- "integrated"

all_samples.integrated <- ScaleData(all_samples.integrated, verbose = T)

all_samples.integrated <- RunPCA(all_samples.integrated, npcs = 50, verbose = T)

# Visualize Variable features
top15 <- head(VariableFeatures(all_samples.integrated), 15)

top15

# Determine the dimentionality

ElbowPlot(all_samples.integrated, ndims = 50)

data <- all_samples.integrated

data <- RunHarmony(data, "HTO_classification")

DimPlot(data, reduction = "harmony", pt.size = 1, group.by = "HTO_classification")

DimPlot(data, reduction = "harmony", pt.size = 1, group.by = "HTO_classification", split.by = "HTO_classification")

data <- data %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.4) %>% identity()

g1 <- DimPlot(data, reduction = "umap", pt.size = 1, split.by = "HTO_classification")
g2 <- DimPlot(data, reduction = "umap", pt.size = 0.5, split.by = "orig.ident")
g3 <- DimPlot(data, reduction = "umap", pt.size = 0.75)
g4 <- DimPlot(data, reduction = "umap", pt.size = 0.75, label = T)
g5 <- DimPlot(data, reduction = "umap", pt.size = 0.5, split.by = "orig.ident", label = T)

DimPlot(data, reduction = "umap", pt.size = 0.5, group.by = "HTO_classification")

ggsave(plot = g2, filename = "Umap split treatment res0.4 dims40.png", path = ".../samples_15_13_1_5/", width = 6.5, height = 5, dpi = 300)

ggsave(plot = g1, filename = "Umap split sample res0.4 dims40.png", path = ".../samples_15_13_1_5/", width = 15, height = 5, dpi = 300)

ggsave(plot = g3, filename = "Umap res0.4 dims40.png", path = ".../samples_15_13_1_5/", width = 6.5, height = 5, dpi = 300)

ggsave(plot = g4, filename = "Umap lable res0.4 dims40.png", path = ".../samples_15_13_1_5/", width = 6.5, height = 5, dpi = 300)

ggsave(plot = g5, filename = "Umap lable split res0.4 dims40.png", path = ".../samples_15_13_1_5/", width = 6.5, height = 5, dpi = 300)

saveRDS(data, file = ".../samples_15_13_1_5/cNF_harmony_integrated.RDS")