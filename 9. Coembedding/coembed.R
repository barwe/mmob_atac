library(ggplot2)
library(cowplot)
library(stringr)

{
    d = readRDS("E:/project/single_cell/mmob_atac/01 AlignmentQC/compared_with_scRNA/medi_coembed_cca.rds")
    d = SplitObject(d, split.by = "tech")
    rna = d$rna
    rm(d)
    rna$orig.ident[is.na(rna$orig.ident)] = "NA"
    rna = SplitObject(rna, split.by = "orig.ident")$WT1
    DefaultAssay(rna)
    rna[['ATAC']] = NULL
    rna[['ACTIVITY']] = NULL
    meta = rna@meta.data
    meta = meta[,!is.na(meta[1,])]
    meta = meta[,!str_ends(colnames(meta), "ATAC")]
    meta = meta[,!str_ends(colnames(meta), "ACTIVITY")]
    meta$res.1.6 = NULL
    ct = meta$ClusterName
    kws = c("OEC", "Astro", "N", "OPC", "MicroG")
    unique(ct)
    for (.i in kws) {
        ct[str_starts(ct, .i)] = .i
    }
    ct[!ct %in% kws] = "else"
    meta$celltype2 = ct
    rna@meta.data = meta
    saveRDS(meta, "metaOfRNA.rds")
    saveRDS(rna, "scRNAseq_WT1_seurat.rds")
} # make rna

atac = readRDS("scATACseq_10x_seurat.rds")
rna = readRDS("scRNAseq_WT1_seurat.rds")

cts = c("N", "OEC", "OPC", "Astro", "MicroG")

intersective.genes = intersect(rownames(atac), rownames(rna))
#atac = atac[intersective.genes, atac$celltype %in% names(table(atac$celltype))[table(atac$celltype) >= 0.008 * ncol(atac)]]
atac = atac[intersective.genes,]
#rna = rna[intersective.genes, rna$celltype %in% names(table(rna$celltype))[table(rna$celltype) >= 0.008 * ncol(rna)]]
rna = rna[intersective.genes,]

atac = NormalizeData(atac)
atac = FindVariableFeatures(atac, nfeatures = 3000)
rna = NormalizeData(rna)
rna = FindVariableFeatures(rna, nfeatures = 3000)

merged = merge(atac, rna)
merged = NormalizeData(merged)
merged = ScaleData(merged, features = VariableFeatures(rna))

merged.ls = SplitObject(merged, split.by = "tech")
atac = merged.ls$atac
rna = merged.ls$rna
rm(merged, merged.ls)


#### Mean

atac.ls = SplitObject(atac, split.by = "celltype2")
atac.ls = lapply(atac.ls, GetAssayData, assay = "ACTIVITY")
atac.ls = lapply(atac.ls, as.matrix)
atac.ls = sapply(atac.ls, rowMeans)

rna.ls = SplitObject(rna, split.by = "celltype2")
rna.ls = lapply(rna.ls, GetAssayData, assay = "RNA")
rna.ls = lapply(rna.ls, as.matrix)
rna.ls = sapply(rna.ls, rowMeans)

library(corrplot)
col = colorRampPalette(c("red","white", "blue"))

spm = cor(atac.ls, rna.ls, method = "spearman")
spm = spm[cts, cts]
saveRDS(spm, "Spearman_coef_mtx.rds")
pdf("Spearman_coef.pdf", width = 12, height = 6)
corrplot(
    spm, 
    method = 'pie',
    type = "full",
    col = col(10),
    cl.pos = "n", #颜色图例不展示
    addCoef.col = rgb(0.2,0.2,0.2), #添加相关系数并设定颜色
    tl.col = "black",
    is.corr = F,
    cl.lim = c(0,1)
)
dev.off()



##################################### Co-Embedding #############################################

library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)

atac = readRDS("scATACseq_10x_seurat.rds")
rna = readRDS("scRNAseq_WT1_seurat.rds")

atac = FindVariableFeatures(atac, nfeatures = 3000)
rna = FindVariableFeatures(rna, nfeatures = 3000)
#saveRDS(rna, "scRNAseq_WT1_seurat.rds")

transfer.anchors = FindTransferAnchors( # slow
    reference = rna, 
    query = atac, 
    features = VariableFeatures(rna), 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    #reduction = "pcaproject"
    reduction = "cca"
)

celltype.predictions = TransferData(
    anchorset = transfer.anchors, 
    refdata = rna$ClusterName, 
    weight.reduction = "cca",
    l2.norm = T,
    k.weight = 10,
    sd.weight = 3,
    verbose = F
)

saveRDS(celltype.predictions, file = "celltype.predictions.rds")
saveRDS(transfer.anchors, file = "transfer.anchors.rds")


atac = readRDS("scATACseq_10x_seurat.rds")
celltype.predictions = readRDS("celltype.predictions.rds")
atac = AddMetaData(atac, metadata = celltype.predictions)
table(celltype.predictions$prediction.score.max > 0.5)


pdf("prediction.score.max.pdf", width = 10, height = 6)
hist(atac$prediction.score.max, main = "Prediction Scores", xlab = "max prediction score")
abline(v = 0.5, col = "red")
dev.off()

rna = readRDS("scRNAseq_WT1_seurat.rds")
atac = readRDS("scATACseq_10x_seurat.rds")
celltype.predictions = readRDS("celltype.predictions.rds")
atac = AddMetaData(atac, metadata = celltype.predictions)
#atac = subset(atac, subset = prediction.score.max > 0.5); print(dim(atac))
atac$predicted.id = factor(atac$predicted.id, levels = unique(rna$ClusterName)) # to make the colors match

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]
transfer.anchors = readRDS("transfer.anchors.rds")
imputation = TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = "cca")
saveRDS(imputation, file = "imputation.rds")

imputation = readRDS("imputation.rds")
atac[["RNA"]] = imputation

rna$celltype = rna$ClusterName
atac$ClusterName = atac$celltype
atac$celltype2 = atac$celltype

coembed = merge(x = rna, y = atac)
coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed = RunUMAP(coembed, dims = 1:30)

#coembed$celltype3 = ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p = DimPlot(coembed, group.by = "tech")
ggsave("coembed_tech.pdf", p, width = 10, height = 10, device = "pdf")
ggsave("coembed_tech.jpg", p, width = 10, height = 10, device = "jpg")

saveRDS(coembed, file = "coembed.rds")
coembed = readRDS("coembed.rds")

coembed$cts = as.factor(paste0(coembed$tech, '_', coembed$celltype2))
p = DimPlot(coembed, group.by = "cts")
ggsave("coembed_celltype_le.pdf", p, width = 10, height = 10, device = "pdf")
ggsave("coembed_celltype_le.jpg", p, width = 10, height = 10, device = "jpg")
p = DimPlot(coembed, group.by = "cts",label = T)+NoLegend()
ggsave("coembed_celltype_nole.pdf", p, width = 10, height = 10, device = "pdf")
ggsave("coembed_celltype_nole.jpg", p, width = 10, height = 10, device = "jpg")

coembed.ls = SplitObject(coembed, split.by = "tech")
atac = coembed.ls$atac
rna = coembed.ls$rna
rna = rna[,rna$celltype2 != "else"]

atac$celltype_ = factor(atac$celltype2, levels = c('N', 'OEC', 'OPC', 'Astro', 'MicroG'))
rna$celltype_ = factor(rna$celltype2, levels = c('N', 'OEC', 'OPC', 'Astro', 'MicroG'))


p11 = DimPlot(atac, group.by = "celltype_", label = T, repel = T) + NoLegend() + ggtitle("scATAC-seq")
p21 = DimPlot(rna, group.by = "celltype_", label = T, repel = T) + NoLegend() + ggtitle("scRNA-seq")
ggsave("sep_atac_nole.pdf", p11, width = 6, height = 6, device = "pdf")
ggsave("sep_rna_nole.pdf", p21, width = 6, height = 6, device = "pdf")
ggsave("sep_atac_nole.jpg", p11, width = 6, height = 6, device = "jpg")
ggsave("sep_rna_nole.jpg", p21, width = 6, height = 6, device = "jpg")

p12 = DimPlot(atac, group.by = "celltype_", label = F, repel = T) + ggtitle("scATAC-seq")
p22 = DimPlot(rna, group.by = "celltype_", label = F, repel = T) + ggtitle("scRNA-seq")
ggsave("sep_atac_le.pdf", p12, width = 8, height = 6, device = "pdf")
ggsave("sep_rna_le.pdf", p22, width = 8, height = 6, device = "pdf")
ggsave("sep_atac_le.jpg", p12, width = 8, height = 6, device = "jpg")
ggsave("sep_rna_le.jpg", p22, width = 8, height = 6, device = "jpg")


predicted.id = atac$predicted.id
kws = c("OEC", "Astro", "N", "OPC", "MicroG")
unique(predicted.id)
for (.i in kws) {
    predicted.id[str_starts(predicted.id, .i)] = .i
}
predicted.id[!predicted.id %in% kws] = "else"


tb = as.data.frame(table(atac$celltype_, predicted.id))
colnames(tb) = c("raw", "pred", "freq")
tb = reshape2::dcast(tb, raw~pred, value.var = "freq")
tb = `rownames<-`(tb[,-1], tb[,1])
tb = tb/rowSums(tb)
orders = c("N", "OEC", "OPC", "Astro", "MicroG")
#tb = tb[orders, c(orders, 'else')]
tb = tb[orders, orders]

tb2 = as.matrix(tb)
col = colorRampPalette(c("red","white", "blue"))
col = colorRampPalette(c("red","white", rgb(54/255, 91/255, 168/255)))


pdf("cell_counts.pdf", width = 8, height = 8)
corrplot::corrplot(
    tb2, 
    method = 'color',
    type = "full",
    col = col(10),
    cl.pos = "n", #颜色图例不展示
    addCoef.col = rgb(0,0,0), #添加相关系数并设定颜色
    tl.col = "black",
    is.corr = F,
    cl.lim = c(0,1)
)
dev.off()
