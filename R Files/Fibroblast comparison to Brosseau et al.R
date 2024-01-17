data <- harmony_all_samples_integrated_cell_type

data <- subset(data, idents = "Fibroblasts")

DefaultAssay(data) <- "RNA"
FeaturePlot(data, features = "FN1")

VlnPlot(data, features = "FN1", split.by = "orig.ident")
VlnPlot(data, features = "COL1A1", split.by = "orig.ident")
VlnPlot(data, features = "COL1A2", split.by = "orig.ident")

VlnPlot(data, features = "CCL2", split.by = "orig.ident")

Fibroblast <- c("CXCL3", "LAMA4", "IL6", "COMP", "HSPA6", "CXCL2", "C11ORF96", "CCL2", "HSPH1", "SFRP2", "PPAP2B", "MT-CO1", "FBLN1", "MT-CO3", "CD9", "MT-ND1", "MT-CYB", "MFAP4", "CTHRC1", "SEPT7", "CEBPD", "H3F3A", "DNAJB1", "HSPA1A", "HSPA1B", "APOD", "RPL41", "RSP17", "CXCL8", "RPL39" )

Idents(data) <- "orig.ident"

g2 <- DotPlot(data, assay = "RNA", features = Fibroblast, scale = T ) + theme_classic() + theme(axis.text = element_text(face = "bold")) + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(face = "bold", size = 10, angle = 0, colour = "black"))+ theme(axis.text.x = element_text(hjust = 1, vjust = 0.5)) + scale_colour_gradient2(low = "blue", mid = "white", high = "red") + scale_size(range = c(1, 6)) + xlab("Gene") + ylab("Genotype") + ggtitle("Brosseau et al comp fibroblast") +theme(title = element_text(face = "bold", size = 12, angle = 0, colour = "black")) + theme(legend.text = element_text(family = 'Helvetica', size = 10)) + theme(legend.title = element_text(family = 'Helvetica', size = 10)) + theme(legend.key.size = unit(0.3, 'cm')) + theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = g2, filename = "Brosseau et al comp fibroblast.png", path = ".../Fibroblast comparion Brosseau et al/", width = 7, height = 3, dpi = 300)