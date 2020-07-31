library(GenomicRanges)
library(Mus.musculus)
library(SnapATAC)
library(Seurat)
library(ggplot2)

# 从小鼠基因组数据注释包Mus.musculus中提取所有的promoter信息
# TSS上游2000和下游500的区域被认为是promoter区域
all_genes = genes(Mus.musculus, "SYMBOL")
all_gene_TSS = resize(all_genes, 1)
promoter = promoters(all_gene_TSS, 2000, 500)
promoter$name = promoter$SYMBOL

# 从scATAC数据集中计算每一个promtoer区域上的reads计数，因为使用的bmat
# 如果使用pmat，就是对peak进行计数，可以考虑一下……
# 计数结果保存到 promoter_mat.rds 文件中，是一个 numCells x numPromoters 的矩阵
xsp = readRDS("D:/project/single_cell/mmob_atac/xsp.rds")
xsp = createGmatFromMat(xsp, "bmat", promoter, do.par = TRUE)
saveRDS(xsp@gmat, "promoter_mat.rds")
xsp = createGmatFromMat(xsp, "pmat", promoter, do.par = TRUE)

# 从 promoter_mat.rds 中筛选出注释为OEC的所有细胞
# 需要导入每个细胞的注释信息。结果保存在 OEC_promote_mat 中
OEC_promote_mat = local({
    promote_mat = readRDS("promoter_mat.rds")
    xsp_cluster = readRDS("D:/project/single_cell/mmob_atac/xsp__cluster.rds")
    promote_mat[which(xsp_cluster == 10),]
})

## psum: accessible degree of each promoter
# 用所有细胞中的计数和来表示某一个promoter的开放程度
psum = sort(colSums(OEC_promote_mat), decreasing = T)

# expression varification of pstrong in scRNA
# 在scRNA数据集上验证开放程度较高的gene的表达情况
rna = readRDS("../rna.rds")

## High accessible analysis ## 

dir.create("RNA_VARI_promoters_high_10")
# 这里多取几个gene，因为有些基因在scRNA数据集上没有表达
pstrong = names(psum)[1:50]
# 在scRNA数据集上绘制小提琴图，肉眼观察gene表达情况
.=lapply(seq_along(pstrong), function(idx, start=1){
    if (idx >= start) {
        pr = pstrong[idx]
        f = sprintf("RNA_VARI_promoters_high_10/%s_%s.png", idx, pr)
        if (!file.exists(f) && pr %in% rownames(rna)) {
            p = VlnPlot(rna, pr, pt.size = 0, group.by = "ClusterName")
            ggsave(f, p, device = "png", width = 12, height = 8)
            message(sprintf("Save to %s", f))
        }
    }
})

# top10 promoters whose gene were detected in scRNA
# 手动挑选出在scRNA数据集OEC上有表达信号的前十个gene
# 并不是pstrong里面的每一个gene都能在scRNA数据集上验证表达
# 小提琴图的命名方式为: [pstrong中的序号]_[gene名称].png
# 可以看到验证出的第10个在scRNA数据集上表达的gene是pstrong中的第35个
hight_top10 = sapply(list.files("RNA_VARI_promoters_high_10"),function(s)strsplit(strsplit(s,'\\.')[[1]][1],'_')[[1]][2])

# 转录因子和靶基因数据库使用TRRUST查询
# 对hight_top10中的每一个gene，在trrust数据库中查询它的转录因子
# 过滤掉没有查询到转录因子的gene，结果保存在tfs中
trrust = readRDS("TRRUST.rds")
tfs = lapply(hight_top10, function(gene) trrust$tf[trrust$target == gene])
names(tfs) = hight_top10
tfs = tfs[sapply(tfs, length) != 0]

# 验证每一个gene的tfs在scRNA数据集上的表达情况
.opdir = "RNA_VARI_factors_of_target"
if(!dir.exists(.opdir))dir.create(.opdir)
for(.gene in names(tfs)){
    .tfs = tfs[[.gene]]
    for(.tf in .tfs){
        f = sprintf("%s/%s_%s.png", .opdir, .gene, .tf)
        if (!file.exists(f) && .tf %in% rownames(rna)) {
            p = VlnPlot(rna, .tf, pt.size = 0, group.by = "ClusterName")
            ggsave(f, p, device = "png", width = 12, height = 8)
            message(sprintf("Save to %s", f))
        }
    }
}

# 手动筛选后重载有表达的tf和它的target
.df = sapply(list.files(.opdir),function(s)strsplit(strsplit(s,'\\.')[[1]][1],'_')[[1]])
rownames(.df) = c("target", "tf")
colnames(.df) = NULL
# 从tfs中提取出有tf表达的target
tfs.filter = tfs[.df["target",]]
# 将其转换为一个dataframe
tfs.df = lapply(names(tfs.filter), function(t){
    tfs = unique(tfs.filter[[t]])
    num = length(tfs)
    targets = rep(t, num)
    df = data.frame(target=targets, tf=tfs, stringsAsFactors = F)
    df$detected = rep(0, num)
    df$detected[df$tf %in% .df["tf",][.df["target",] == t]] = 1
    df
})
df = Reduce(rbind, tfs.df)
# 导出csv文件，cytoscape展示
df$tf_id = paste0(df$tf, "->", df$target)
write.csv(df, "figure_1.csv", quote = F, row.names = F)


## Figure 2: low accessible analysis ##
# 这里的计数最小值很集中，无从判断top10
trrust = readRDS("TRRUST.rds")
.opdir = "RNA_VARI_promoters_low2"
if (!dir.exists(.opdir)) dir.create(.opdir)
psum = sort(colSums(OEC_promote_mat))
psum = psum[psum > 0]
.=lapply(names(psum)[300:350], function(gene){
    fp = sprintf("%s/%s.png", .opdir, gene)
    if (!file.exists(fp) && gene %in% rownames(rna)) {
        p = VlnPlot(rna, gene, pt.size = 0, group.by = "ClusterName")
        ggsave(fp, p, device = "png", width = 12, height = 8)
        message(sprintf("Save to %s", fp))
    }
})

low = sapply(list.files(.opdir),function(s)strsplit(s, "\\.")[[1]][1])
tfs = lapply(low, function(gene) trrust$tf[trrust$target == gene])
names(tfs) = low
tfs = tfs[sapply(tfs, length) != 0]
tfs

# 验证每一个gene的tfs在scRNA数据集上的表达情况
.opdir = "RNA_VARI_factors_of_target_low2"
if(!dir.exists(.opdir))dir.create(.opdir)
for(.gene in names(tfs)){
    .tfs = tfs[[.gene]]
    for(.tf in .tfs){
        f = sprintf("%s/%s_%s.png", .opdir, .gene, .tf)
        if (!file.exists(f) && .tf %in% rownames(rna)) {
            p = VlnPlot(rna, .tf, pt.size = 0, group.by = "ClusterName")
            ggsave(f, p, device = "png", width = 12, height = 8)
            message(sprintf("Save to %s", f))
        }
    }
}

.df = sapply(list.files(.opdir),function(s)strsplit(strsplit(s,'\\.')[[1]][1],'_')[[1]])
rownames(.df) = c("target", "tf")
colnames(.df) = NULL
# 从tfs中提取出有tf表达的target
tfs.filter = tfs[.df["target",]]
# 将其转换为一个dataframe
tfs.df = lapply(names(tfs.filter), function(t){
    tfs = unique(tfs.filter[[t]])
    num = length(tfs)
    targets = rep(t, num)
    df = data.frame(target=targets, tf=tfs, stringsAsFactors = F)
    df$detected = rep(0, num)
    df$detected[df$tf %in% .df["tf",][.df["target",] == t]] = 1
    df
})
df = Reduce(rbind, tfs.df)
# 导出csv文件，cytoscape展示
df$tf_id = paste0(df$tf, "->", df$target)
write.csv(df, "figure_2_2.csv", quote = F, row.names = F)





