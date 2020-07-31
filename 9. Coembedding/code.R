library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

#### Preprocess ATAC

atac.peakdata = readRDS("peak_cell_dgCMatrix.rds")
rownames(atac.peakdata) = readRDS("peak-name.rds")
colnames(atac.peakdata) = readRDS("barcode.rds")

activity.matrix = CreateGeneActivityMatrix(
    peak.matrix = atac.peakdata, 
    annotation.file = "Mus_musculus.GRCm38.99.chr.gtf", 
    seq.levels = c(1:19, "X", "Y"), 
    upstream = 2000, verbose = T)

atac = CreateSeuratObject(counts = atac.peakdata, assay = "ATAC", project = "10x_scATAC")
atac[["ACTIVITY"]] = CreateAssayObject(counts = activity.matrix)
atac.meta = readRDS("xsp__meta2.rds")
atac = AddMetaData(atac, metadata = atac.meta)
atac$tech = "atac"

DefaultAssay(atac) <- "ACTIVITY"
atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)

DefaultAssay(atac) <- "ATAC"
var.peaks = names(which(Matrix::rowSums(atac) > 300))
message("peak var/all ratio:", length(var.peaks)/nrow(atac))
VariableFeatures(atac) <- var.peaks
atac <- RunLSI(atac, n = 50, scale.max = NULL) # slow
atac <- RunUMAP(atac, reduction = "lsi", dims = 1:50)

saveRDS(activity.matrix, file = "medi_activity.matrix.rds")
saveRDS(atac, file = "medi_atac_seurat.rds")


#### Preprocess RNA

rna.data = fread("GSE121891_OB_6_runs.raw.dge.csv", quote = F)
rna.data = `rownames<-`(as.matrix(rna.data[,-1]), rna.data[,1]$V1)

rna = CreateSeuratObject(counts = rna.data, assay = "RNA", project = "scRNA")
rna.meta = read.table("GSE121891_OB_metaData_seurat.csv", sep=",", header=T, stringsAsFactors=F)
rna.meta = `rownames<-`(rna.meta[,-1], rna.meta[,1])
rna.meta = rna.meta[colnames(rna),]
rna = AddMetaData(rna, metadata = rna.meta)
rna$celltype = stringr::str_match(rna$ClusterName, "(\\D+).*")[,2]
rna$tech = "rna"

rna <- FindVariableFeatures(rna)
rna <- NormalizeData(rna)
rna <- ScaleData(rna)

rna = RunPCA(rna, features = VariableFeatures(rna))
rna = JackStraw(rna, num.replicate = 100) # slow
{
    pdf("rna-ElbowPlot.pdf")
    ElbowPlot(rna)
    dev.off()
}

saveRDS(rna, file = "medi_rna_seurat.rds")


#### 

atac <- readRDS("medi_atac_seurat.rds")
rna <- readRDS("medi_rna_seurat.rds")
rna <- RunUMAP(rna, dims = 1:20)
{
    p1 <- DimPlot(atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
    p2 <- DimPlot(rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
    pdf("UMAP-2D.pdf", width = 12, height = 6)
    p1 + p2
    dev.off()
} # plot
transfer.anchors <- FindTransferAnchors( # slow
    reference = rna, 
    query = atac, 
    features = VariableFeatures(rna), 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca")

celltype.predictions <- TransferData(
    anchorset = transfer.anchors, 
    refdata = rna$celltype, 
    weight.reduction = "cca",
    l2.norm = T,
    k.weight = 10,
    sd.weight = 3,
    verbose = F
)
atac = AddMetaData(atac, metadata = celltype.predictions)
table(atac$prediction.score.max > 0.5)

saveRDS(rna, file = "medi_rna_seurat.rds")
saveRDS(transfer.anchors, file = "medi_transfer.anchors.rds")
saveRDS(celltype.predictions, file = "medi_celltype.predictions.rds")


#### 

atac = readRDS("medi_atac_seurat.rds")
rna = readRDS("medi_rna_seurat.rds")
celltype.predictions = readRDS("medi_celltype.predictions.rds")
{
    pdf("prediction.score.max.pdf", width = 10, height = 6)
    hist(atac$prediction.score.max)
    abline(v = 0.5, col = "red")
    dev.off()
} # plot

## view the predicted cell types on a UMAP representation of the scATAC-seq data and 
## find that the transferred labels are highly consistent with the UMAP structure.
atac.filtered <- subset(atac, subset = prediction.score.max > 0.5)
atac.filtered$predicted.id <- factor(atac.filtered$predicted.id, levels = unique(rna$celltype))  # to make the colors match
{
    p1 <- DimPlot(atac.filtered, group.by = "predicted.id", label = T, repel = TRUE) + 
        ggtitle("scATAC-seq cells") + 
        NoLegend() + scale_colour_hue(drop = FALSE)
    p2 <- DimPlot(rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
        NoLegend()
    p = cowplot::plot_grid(p1, p2, ncol = 2)
    ggsave("medi_1p.pdf", p, width = 12, height = 6)
} # plot


#### Co-embedding

atac = readRDS("medi_atac_seurat.rds")
rna = readRDS("medi_rna_seurat.rds")
transfer.anchors = readRDS("medi_transfer.anchors.rds")

# version1: weight.reduction = atac[["lsi"]]
{
    genes.use <- VariableFeatures(rna)
    refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]])
    atac[["RNA"]] <- imputation
    coembed <- merge(x = rna, y = atac) # slow
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
    
    saveRDS(coembed, file = "medi_coembed_lsi.rds")
    
    coembed = readRDS("medi_coembed_lsi.rds")
    coembed <- RunUMAP(coembed, dims = 1:30)
    coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
    
    coembed.filter = coembed[,!is.na(coembed$celltype)]
    suffix = ifelse(coembed.filter$tech == "atac", "_scATAC", "_scRNA")
    celltype = coembed.filter$celltype
    celltype[celltype == "Mφ"] = "Mphi"
    coembed.filter$celltype2 = paste0(celltype, suffix)
    
    p1 <- DimPlot(coembed.filter, group.by = "tech")
    p2 <- DimPlot(coembed.filter, group.by = "celltype2", label = T, repel = TRUE)
    ggsave("medi_coembed1_lsi.pdf", p1, width = 12, height = 10)
    ggsave("medi_coembed2_lsi.pdf", p2, width = 12, height = 8)
    #ggsave("medi_coembed_lsi.pdf", cowplot::plot_grid(p1, p2, ncol = 2), width = 12, height = 6)
    
    tmp = coembed.filter[,match(unique(coembed.filter$celltype2), coembed.filter$celltype2)]
    tmp = DimPlot(tmp, group.by = "celltype2", label = T, repel = TRUE)
    ggsave("medi_coembed3_lsi.pdf", tmp, width = 12, height = 8)
    
    saveRDS(coembed, file = "medi_coembed_lsi.rds")
} 

# version2: weight.reduction = "cca"
{
    genes.use <- VariableFeatures(rna)
    refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = "cca")
    atac[["RNA"]] <- imputation
    coembed <- merge(x = rna, y = atac)
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
    
    saveRDS(coembed, file = "medi_coembed_cca.rds")
    
    coembed = readRDS("medi_coembed_cca.rds")
    coembed <- RunUMAP(coembed, dims = 1:30)
    coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
    
    p1 <- DimPlot(coembed, group.by = "tech")
    p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
    ggsave("medi_coembed1_cca.pdf", p1, width = 10, height = 10)
    ggsave("medi_coembed2_cca.pdf", p2, width = 10, height = 10)
    
    ggsave("medi_coembed_cca.pdf", cowplot::plot_grid(p1, p2, ncol = 2), width = 12, height = 6)

    saveRDS(coembed, file = "medi_coembed_cca.rds")
}


coembed.ls = SplitObject(coembed.filter, split.by = "tech")

ps = lapply(coembed.ls, DimPlot, group.by = "celltype2", label = T, repel = T)

ggsave("soembed_split.pdf", plot_grid(plotlist = ps, ncol = 2), width = 12, height = 5)








