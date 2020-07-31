suppressMessages(library(Matrix))
suppressMessages(library(ChIPseeker))

peakAnno = readRDS("peakAnno.ls.rds")
genes = lapply(peakAnno, function(anno) {
    pos = anno@detailGenomicAnnotation
    pos = pos$Exon & pos$Promoter
    id = anno@anno@elementMetadata@listData$geneId[pos]
    symbol = anno@anno@elementMetadata@listData$SYMBOL[pos]
    names(symbol) = as.integer(id)
    symbol = symbol[unique(names(symbol))]
    return(symbol)
})
genes.counts = sapply(genes, length)


tfs = read.table("TF.all.txt", stringsAsFactors = F, sep = '\t', header = T)
tfs.counts = sapply(split(tfs$TF, tfs$Cluster), length)

library(ggplot2)
cinfo = readRDS("cinfo.rds")
cell.type = cinfo$cell.type[order(cinfo$ordered.num)]

wdf = data.frame(
    x = paste0('C', 1:21, '\n', cell.type), 
    tf =  tfs.counts, 
    ta = genes.counts
)
saveRDS(wdf, "wdf.rds")
#ldf = reshape2::melt(wdf, values = c("target", "factor"))


wdf = readRDS("wdf.rds")
df = data.frame(
    order = 1:21,
    x = wdf$x[order(wdf$tf, decreasing = T)],
    y = sort(wdf$tf, decreasing = T)
)

ggplot(df, aes(order, y, fill = "red")) + 
    geom_bar(stat = "identity") +
    labs(x = "", y = "Counts") +
    ggtitle("Hist of enriched factors") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks=1:21, labels = df$x)+
    guides(fill = FALSE) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks.x = element_blank(),
    )
ggsave("TF_counts.pdf", width = 6, height = 4)


df = data.frame(
    order = 1:21,
    x = wdf$x[order(wdf$ta, decreasing = T)],
    y = sort(wdf$ta, decreasing = T)
)
ggplot(df, aes(order, y, fill = "red")) + 
    geom_bar(stat = "identity") +
    labs(x = "", y = "Counts") +
    ggtitle("Hist of candidated targets") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks=1:21, labels = df$x)+
    guides(fill = FALSE) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks.x = element_blank(),
    )

ggsave("target_counts.pdf", width = 6, height = 4)
