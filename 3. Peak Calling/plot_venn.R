## All peaks which can be annotated to their nearest gene
peakAnno.ls <- readRDS("peakAnno.ls.rds")
cinfo <- readRDS("D:/project/single_cell/mmob_atac/cinfo.rds")

## Annotated position status matrix list
pos.status.ls = lapply(peakAnno.ls, function(d) d@detailGenomicAnnotation)

## Annotated gene ids of each cluster
geneId.ls = lapply(peakAnno.ls, function(d) d@anno@elementMetadata@listData$geneId)

## Only annotating within promoter
output_dir = "venn_data_in_promoter"
dir.create(output_dir)
gids = lapply(1:21, function(i) unique(geneId.ls[[i]][pos.status.ls[[i]]$Promoter]))
names(gids) = 1:21
cts = unique(cinfo$cell.type)
ct.ls = lapply(cts, function(s) {
    cluster.ids = cinfo$ordered.num[cinfo$cell.type == s]
    genes = unique(Reduce(c, gids[cluster.ids]))
    df = data.frame(x=genes,stringsAsFactors = F)
    fp = sprintf("%s/%s", output_dir, s)
    write.table(df, fp, quote = F, col.names = F, row.names = F)
})

## Within promoter and gene
output_dir = "venn_data_in_both"
dir.create(output_dir)
gids = lapply(1:21, function(i) unique(geneId.ls[[i]][pos.status.ls[[i]]$Promoter & pos.status.ls[[i]]$genic]))
names(gids) = 1:21
cts = unique(cinfo$cell.type)
ct.ls = lapply(cts, function(s) {
    cluster.ids = cinfo$ordered.num[cinfo$cell.type == s]
    genes = unique(Reduce(c, gids[cluster.ids]))
    df = data.frame(x=genes,stringsAsFactors = F)
    fp = sprintf("%s/%s", output_dir, s)
    write.table(df, fp, quote = F, col.names = F, row.names = F)
})

## Within gene
output_dir = "venn_data_in_gene"
dir.create(output_dir)
gids = lapply(1:21, function(i) unique(geneId.ls[[i]][pos.status.ls[[i]]$genic]))
names(gids) = 1:21
cts = unique(cinfo$cell.type)
ct.ls = lapply(cts, function(s) {
    cluster.ids = cinfo$ordered.num[cinfo$cell.type == s]
    genes = unique(Reduce(c, gids[cluster.ids]))
    df = data.frame(x=genes,stringsAsFactors = F)
    fp = sprintf("%s/%s", output_dir, s)
    write.table(df, fp, quote = F, col.names = F, row.names = F)
})
