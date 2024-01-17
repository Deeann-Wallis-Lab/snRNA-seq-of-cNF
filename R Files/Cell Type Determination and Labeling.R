DimPlot(cNF_harmony_integrated, pt.size = 0.5, split.by = "orig.ident")

DefaultAssay(cNF_harmony_integrated) <- "RNA"

DimPlot(object = cNF_harmony_integrated, pt.size = 0.5, reduction = "umap", label = T)


# Schwann cells
Schwann_cells <- cNF_Genes_List_for_Cell_Population_Identification$`Schwann Cells`

g <- FeaturePlot(cNF_harmony_integrated, features = Schwann_cells)

ggsave(plot = g,
       filename = "Schwann_cells.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)


#Endothelial cells
Endothelial_cells <- cNF_Genes_List_for_Cell_Population_Identification$`Endothelial Cells`

g <- FeaturePlot(cNF_harmony_integrated, features = Endothelial_cells)

ggsave(plot = g,
       filename = "Endothelial_cells.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

#Epithelial cells
Epithelial_cells <- cNF_Genes_List_for_Cell_Population_Identification$`Epithelial Cells`

g <- FeaturePlot(cNF_harmony_integrated, features = Epithelial_cells)

ggsave(plot = g,
       filename = "Epithelial_cells.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

#Fibroblasts cells
Fibroblasts <- cNF_Genes_List_for_Cell_Population_Identification$Fibroblasts

g <- FeaturePlot(cNF_harmony_integrated, features = Fibroblasts)

ggsave(plot = g,
       filename = "Fibroblasts.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

#Pericytes cells
Pericytes <- cNF_Genes_List_for_Cell_Population_Identification$Pericytes

g <- FeaturePlot(cNF_harmony_integrated, features = Pericytes)

ggsave(plot = g,
       filename = "Pericytes.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

#Macrophages
Macrophages <- cNF_Genes_List_for_Cell_Population_Identification$Macrophages

Macrophage_features <- c("MRC1", "CD163", "ITGAM", "RBPJ", "MS4A4A", "MS4A6A", "CD16", "C1Qa", "MSR1","FCER1G")

g <- FeaturePlot(cNF_harmony_integrated, features = Macrophages)

g1 <- FeaturePlot(cNF_harmony_integrated, features = Macrophage_features)

ggsave(plot = g,
       filename = "Macrophages2.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

ggsave(plot = g1,
       filename = "Macrophages.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

#Mast cells
Mast_cells_1 <- cNF_Genes_List_for_Cell_Population_Identification$`Mast Cells`

Mast_Cell_2 <- c('KIT', "LTC4S", "TPSAB1", "IL1RL1", "TPSB2", " CMA1", "FCER2", "FCGR2A", "CD33", "PTPRC")

g <- FeaturePlot(cNF_harmony_integrated, features = Mast_Cell_2)

g1 <- FeaturePlot(cNF_harmony_integrated, features = Mast_cells_1)

ggsave(plot = g,
       filename = "Mast_cells2.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

ggsave(plot = g1,
       filename = "Mast_cells2.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

# Keratinocytes
Keratinocytes <- cNF_Genes_List_for_Cell_Population_Identification$Keratinocytes

g <- FeaturePlot(cNF_harmony_integrated, features = Keratinocytes)

g1 <- FeaturePlot(cNF_harmony_integrated, features = Keratinocytes)

ggsave(plot = g,
       filename = "Keratinocytes.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

# Epithelial cells
Epithelial_cells <- cNF_Genes_List_for_Cell_Population_Identification$`Epithelial cells`

g <- FeaturePlot(cNF_harmony_integrated, features = Epithelial_cells)

g1 <- FeaturePlot(cNF_harmony_integrated, features = Epithelial_cells)

ggsave(plot = g,
       filename = "Epithelial_cells_2.png",
       path = ".../cell type/",
       width = 15,
       height = 15,
       dpi = 300)

# Melanocytes
Melanocytes <- cNF_Genes_List_for_Cell_Population_Identification$Melanocytes

g <- FeaturePlot(cNF_harmony_integrated, features = Melanocytes)

g1 <- FeaturePlot(cNF_harmony_integrated, features = Melanocytes)

ggsave(plot = g,
       filename = "Melanocytes.png",
       path = ".../samples_15_13_1_5/cell type/",
       width = 15,
       height = 15,
       dpi = 300)


# label cell types in Seurat object

all_samples.integrated <- cNF_harmony_integrated

DimPlot(all_samples.integrated, label = T)

new.cluster.ids <- c("3",	"11",	"2",	"13",	"1",	"10",	"12",	"14",	"8",	"6",	"4",	"9",	"7",	"5")

names(new.cluster.ids) <- levels(all_samples.integrated)

all_samples.integrated <- RenameIdents(all_samples.integrated, new.cluster.ids)

levels(x=all_samples.integrated)<- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

DimPlot(all_samples.integrated, label = T)

new.cluster.ids <- c("Endothelial cells 1",	"Endothelial cells 1",	"Endothelial cells 1",	"Endothelial cells 1",	"Endothelial cells 1",	"Keratinocytes",	"Melanocytes",	"Pericytes",	"Myeloid cells", 	"Fibroblasts",	"Fibroblasts",	"Fibroblasts",	"Schwann cells",	"Endothelial cells 2")


names(new.cluster.ids) <- levels(all_samples.integrated)

all_samples.integrated <- RenameIdents(all_samples.integrated, new.cluster.ids)

g <- DimPlot(all_samples.integrated, label = F)

g1 <- DimPlot(all_samples.integrated, label = F, split.by = "orig.ident")

ggsave(plot = g,
       filename = "cell type.png",
       path = ".../samples_15_13_1_5/",
       width = 6.5,
       height = 5,
       dpi = 300)

ggsave(plot = g1,
       filename = "cell type split .png",
       path = ".../samples_15_13_1_5/",
       width = 8,
       height = 5,
       dpi = 300)


my_cols <- c("Endothelial cells 1" = "#663c8f", "Keratinocytes" = "#9cb241", "Melanocytes" = "#7080d9", "Pericytes" = "#c99135", "Myeloid cells" = "#c868b1", "Fibroblasts" = "#5bbc72", "Schwann cells" =  "#b54361", "Endothelial cells 2" = "#43c9b0")

g <- DimPlot(all_samples.integrated, label = F, cols = my_cols)

g1 <- DimPlot(all_samples.integrated, label = F, split.by = "orig.ident", cols = my_cols)

ggsave(plot = g,
       filename = "cell type recolor.png",
       path = ".../samples_15_13_1_5/",
       width = 8,
       height = 5,
       dpi = 300)

ggsave(plot = g1,
       filename = "cell type split recolor.png",
       path = ".../samples_15_13_1_5/",
       width = 8,
       height = 5,
       dpi = 300)

all_samples.integrated@meta.data[["cell type"]] <- Idents(all_samples.integrated)

saveRDS(all_samples.integrated, file = ".../harmony_all_samples_integrated_cell_type.rds")


n_cells <- FetchData(all_samples.integrated, 
                     vars = c("ident", "HTO_classification")) %>%
  dplyr::count(ident, HTO_classification) %>%
  tidyr::spread(ident, n)

write_xlsx(n_cells,path = ".../harmony cell count by cell type.xlsx")

# calculated frequency of each cell type then imported the data back into R for graphing

bg <- harmony_cell_count_by_cell_type

test2<- melt(bg, id.vars = "group")

g <- ggplot(test2, aes(fill=variable, y=value, x=group)) + 
  geom_bar(position="stack", stat="identity",width = 0.7) + 
  theme_classic() + theme(axis.text = element_text(face = "bold")) + 
  coord_flip() + 
  theme(axis.text = element_text(face = "bold", size = 12, angle = 0, colour = "black"), axis.text.x = element_text(hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 0.5))


ggsave(plot = g,
       filename = "cell type frequency.png",
       path = ".../samples_15_13_1_5/",
       width = 6,
       height = 2,
       dpi = 300)


write_xlsx(n_cells,path = ".../harmony cell count by cell type all_2 endo population.xlsx")



