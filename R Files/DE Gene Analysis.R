# DE genes on subsets 4 samples
# import RDS from Cell Type Determination and Labeling script 

data <- harmony_all_samples_integrated_cell_type

Idents(data) <- "cell type"

Endothelial <- subset(data, idents = "Endothelial cells 1")

Fibroblast <- subset(data, idents = "Fibroblasts")

Schwann_cells <- subset(data, idents = "Schwann cells")

Epithelial_cells_Keratinocytes <- subset(data, idents = "Keratinocytes")

Myeloid_cells <- subset(data, idents = "Myeloid cells")

Pericytes <- subset(data, idents = "Pericytes")

Melanocytes <- subset(data, idents = "Melanocytes")

Endothelial_cells_2 <- subset(data, idents = "Endothelial cells 2")


Idents(data) <- "orig.ident"

all_samples.integrated_DE_treatment_control <- FindMarkers(data, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

all_samples.integrated_DE_treatment_control$genes <- row.names(all_samples.integrated_DE_treatment_control)

#####################################
Idents(Endothelial) <- "orig.ident"

Endothelial_DE_treatment_control <- FindMarkers(Endothelial, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Endothelial_DE_treatment_control$genes <- row.names(Endothelial_DE_treatment_control)

#####################################
Idents(Fibroblast) <- "orig.ident"

Fibroblast_DE_treatment_control <- FindMarkers(Fibroblast, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Fibroblast_DE_treatment_control$genes <- row.names(Fibroblast_DE_treatment_control)

#####################################
Idents(Schwann_cells) <- "orig.ident"

Schwann_cells_DE_treatment_control <- FindMarkers(Schwann_cells, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Schwann_cells_DE_treatment_control$genes <- row.names(Schwann_cells_DE_treatment_control)

#####################################
Idents(Epithelial_cells_Keratinocytes) <- "orig.ident"

Epithelial_cells_Keratinocytes_DE_treatment_control <- FindMarkers(Epithelial_cells_Keratinocytes, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Epithelial_cells_Keratinocytes_DE_treatment_control$genes <- row.names(Epithelial_cells_Keratinocytes_DE_treatment_control)

#####################################
Idents(Myeloid_cells) <- "orig.ident"

Myeloid_cells_DE_treatment_control <- FindMarkers(Myeloid_cells, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Myeloid_cells_DE_treatment_control$genes <- row.names(Myeloid_cells_DE_treatment_control)

#####################################
Idents(Pericytes) <- "orig.ident"

Pericytes_DE_treatment_control <- FindMarkers(Pericytes, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Pericytes_DE_treatment_control$genes <- row.names(Pericytes_DE_treatment_control)

#####################################
Idents(Melanocytes) <- "orig.ident"

Melanocytes_DE_treatment_control <- FindMarkers(Melanocytes, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Melanocytes_DE_treatment_control$genes <- row.names(Melanocytes_DE_treatment_control)


#####################################
Idents(Endothelial_cells_2) <- "orig.ident"

Endothelial_cells_2_DE_treatment_control <- FindMarkers(Endothelial_cells_2, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "treated", ident.2 = "control")

Endothelial_cells_2_DE_treatment_control$genes <- row.names(Endothelial_cells_2_DE_treatment_control)


sheets <- list("all cells" = all_samples.integrated_DE_treatment_control,"Endothelial" = Endothelial_DE_treatment_control, "Fibroblast" = Fibroblast_DE_treatment_control ,  'Schwann_cells' = Schwann_cells_DE_treatment_control,'Epithelial_Keratinocytes' = Epithelial_cells_Keratinocytes_DE_treatment_control, "Myeloid_cells" = Myeloid_cells_DE_treatment_control, "Pericytes"= Pericytes_DE_treatment_control, "Melanocytes" = Melanocytes_DE_treatment_control, "Endothelial_cells_2" = Endothelial_cells_2_DE_treatment_control)

write.xlsx(sheets, file = "...DE_genes_treatmen_control_cell_type.xlsx")





