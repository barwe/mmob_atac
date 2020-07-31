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
.=lapply(1:21, function(i){
    gs = genes[[i]]
    df = data.frame(x = gs)
    fp = sprintf("genes_anno/Cluster_%s.txt", i)
    write.table(df, fp, quote = F, row.names = F, col.names = F)
})

genes.counts = sapply(genes, length)


tfs = read.table("TF.all.txt", stringsAsFactors = F, sep = '\t', header = T)
tfs.counts = sapply(split(tfs$TF, tfs$Cluster), length)

library(ggplot2)
cinfo = readRDS("cinfo.rds")
cell.type = cinfo$cell.type[order(cinfo$ordered.num)]

idx = order(genes.counts, decreasing = T)

wdf = data.frame(
    x = 1:21,
    x.label = paste0('C', 1:21, '\n', cell.type)[idx], 
    #tf =  tfs.counts, 
    ta = genes.counts[idx]
)
#ldf = reshape2::melt(wdf, values = c("target", "factor"))

tfp = ggplot(wdf, aes(x, tf, fill = "red")) + 
    geom_bar(stat = "identity") +
    labs(x = "", y = "Counts") +
    ggtitle("Candidate genes associated with accessibility") +
    #geom_text(aes(label=tf,y=tf+0.4),position=position_dodge(0.9),vjust=0,size=3,alpha=0.5)+
    scale_y_continuous(expand = c(0,0)) +
    #guides(fill = guide_legend(title = "Types")) +
    guides(fill = FALSE) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.position = 'top'
    )


tap = ggplot(wdf, aes(x, ta, fill = "red")) + 
    geom_bar(stat = "identity") +
    labs(x = "", y = "Counts") +
    ggtitle("Hist of accessible genes") +
    #geom_text(aes(label=tf,y=tf+0.4),position=position_dodge(0.9),vjust=0,size=3,alpha=0.5)+
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = 1:21, labels = wdf$x.label)+
    #guides(fill = guide_legend(title = "Types")) +
    guides(fill = FALSE) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.position = 'top'
    )

ggsave("TF_counts.pdf", tfp, width = 6, height = 4)
ggsave("target_counts.pdf", tap, width = 6, height = 4)
