# 需要：
# 1. 每个细胞每个promoter的打分矩阵
# 2. 每个细胞所属的cluster
# 3. 每个cluster的peak注释信息
# 方法：
# 首先挑选出MG的所有细胞的promoter的子矩阵，计算每个promoter的平均打分
# 根据MG的每个peak的注释信息挑选出归属于Promoter或者Exon的peak对应的gene
# 对这些基因的promoter平均打分排序

suppressMessages(library(Matrix))
suppressMessages(library(ChIPseeker))

peakAnno.ls = readRDS("peakAnno.ls.rds")
prmat <- readRDS("promoter_mat.rds")
xsp.cluster <- readRDS("xsp__cluster.rds")

for (i in 1:21) {
    message(i)
    promoters = colMeans(prmat[xsp.cluster == i,])
    peakAnno = peakAnno.ls[[i]]
    pos = peakAnno@detailGenomicAnnotation
    pos = pos$Exon & pos$Promoter
    symbols = peakAnno@anno@elementMetadata@listData$SYMBOL[pos]
    promoters.sel = sort(promoters[symbols], decreasing = T)
    promoters.sel = promoters.sel[unique(names(promoters.sel))]
    promoters.sel = promoters.sel[promoters.sel > 0]
    df = data.frame(
        symbol = names(promoters.sel),
        signal = promoters.sel,
        stringsAsFactors = F
    )
    fp = sprintf("promoterSignal_%s.csv", i)
    write.table(df, fp, sep = '\t', quote = F, row.names = F)
}

