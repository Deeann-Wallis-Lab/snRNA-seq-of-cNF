# CellChat treatment and control comparison analysis Trun mean 0.1 min cells 10

# import Seurat RDS files named harmony_all_samples_integrated_cell_type.RDS

Human_genes_data_slot <- harmony_all_samples_integrated_cell_type

#### treament

options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot)
Idents(Human_genes_data_slot) <- "orig.ident"
data <- subset(Human_genes_data_slot, idents = "treated")
Idents(data) <- "cell type"
DimPlot(data)

cellchat <- createCellChat(object = data, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



# for circle plots

pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

# for heatmaps


pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

#saveRDS(cellchat, file = ".../treated.rds")

#### control

options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot)
Idents(Human_genes_data_slot) <- "orig.ident"
data <- subset(Human_genes_data_slot, idents = "control")
Idents(data) <- "cell type"
DimPlot(data)

cellchat <- createCellChat(object = data, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# for circle plots

pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

# for heatmaps


pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

saveRDS(cellchat, file = ".../control.rds")


# comparison of all Treatment to control

# Define the cell labels to lift up
group.new = levels(treated@idents)
TC <- liftCellChat(control, group.new)


object.list <- list(control = control, treated = treated)

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

TC@netP$pathways

pathways_treated <- treated@netP$pathways
pathways_treated <- as.data.frame(pathways_treated)
colnames(pathways_treated)[1] <- "pathways"

pathways_control <- control@netP$pathways
pathways_control <- as.data.frame(pathways_control)
colnames(pathways_control)[1] <- "pathways"

pathways <- pathways_treated %>% filter(pathways %in% pathways_control$pathways)



for (i in 1:nrow(pathways)){
  
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)  
  
  pathways.show <- pathways[i,]
  
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
  
  vertex.receiver = seq(1,10) 
  
  par(mfrow = c(1,2), xpd=TRUE)
  
  
  
  (for (a in 1:length(object.list)) {
    
    test2 <- netVisual_aggregate(object.list[[a]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[a]))
    
  })
  
  
  test2
  dev.off()
}
dev.off()


object.list <- list(treated = treated, control = control)

dev.off()


pdf(file ="cellchat.pdf", width = 8, height =8)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# incoming and outgoing heatmaps
ct1 <- netAnalysis_signalingRole_heatmap(control, pattern = "outgoing", width = 4,height = 20)
ct2 <- netAnalysis_signalingRole_heatmap(control, pattern = "incoming", width = 4,height = 20)


tt1 <- netAnalysis_signalingRole_heatmap(treated, pattern = "outgoing", width = 4,height = 20)
tt2 <- netAnalysis_signalingRole_heatmap(treated, pattern = "incoming", width = 4,height = 20)



png(filename = paste("outgoing control.png"), width = 4, height =10, units = "in", res = 300)
ct1
dev.off()


png(filename = paste("incoming control.png"), width = 4, height =10, units = "in", res = 300)
ct2
dev.off()


png(filename = paste("outgoing treated.png"), width = 4, height =10, units = "in", res = 300)
tt1
dev.off()


png(filename = paste("incoming treated.png"), width = 4, height =10, units = "in", res = 300)
tt2
dev.off()


