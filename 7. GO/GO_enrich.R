## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"

library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

## Run ChIPseeker::annotatePeak for each cluster
peakAnno.ls = readRDS(paste0(PROJECT_DIR, "/05 peak calling/peakAnno.ls.rds"))

## Extract annotational information into table
peakAnnoTbl.ls = lapply(peakAnno.ls, function(obj) data.frame(obj@anno@elementMetadata@listData, stringsAsFactors = F))
saveRDS(peakAnnoTbl.ls, file = "peakAnnoTbl.ls.rds")

## Extract gene symbols from annotational table
peakAnnoGene.ls = lapply(peakAnnoTbl.ls, function(df) unique(df$SYMBOL))
saveRDS(peakAnnoGene.ls, file = "peakAnnoGene.ls.rds")

## Convert gene symbols into gene entrez id
peakAnnoEntrezid.ls = lapply(peakAnnoGene.ls, function(gs){
    df = bitr(gs, fromType = "SYMBOL", toType = "ENTREZID", 
              OrgDb = org.Mm.eg.db)
    as.vector(df$ENTREZID, "character")
})
saveRDS(peakAnnoEntrezid.ls, file = "peakAnnoEntrezid.ls.rds")

## Run clusterProfiler::enrichGO
run_enrichGO = function(entrezid.ls, num.cores=1){
    if(num.cores==1){
        go.res = lapply(entrezid.ls, enrichGO, OrgDb = org.Mm.eg.db,
                        ont = 'BP', pAdjustMethod = 'BH', readable = TRUE,
                        pvalueCutoff = 1, qvalueCutoff = 1)
    }
    else{
        library(parallel)
        if(num.cores > detectCores()-2)
            stop(sprintf("Current %s available cores and you want to use %s", detectCores(), num.cores))
        mkc = makeCluster(num.cores)
        GO_func = function(entrezids){
            library(org.Mm.eg.db)
            clusterProfiler::enrichGO(entrezids, org.Hs.eg.db, ont = 'BP', 
                                      pAdjustMethod = 'BH', pvalueCutoff = 1, 
                                      qvalueCutoff = 1, readable = TRUE)
        }
        go.res = parLapply(mkc, entrezid.ls, GO_func)
        stopCluster(mkc)
    }
    return(go.res)
}
enrichRes.ls = run_enrichGO(peakAnnoEntrezid.ls, num.cores = 1)
saveRDS(enrichRes.ls, file = "clusterEnrichGO.ls.rds")

## Extract and save tables from result of enrichGO
if(stringr::str_ends(getwd(), work_dir)) setwd("../")

enrichRes.dir = "GO_tables"
if(!dir.exists(enrichRes.dir)) dir.create(enrichRes.dir)

lapply(names(enrichRes.ls), function(s){
    fp = sprintf("%s/GO-enrich_cluster%s.csv",enrichRes.dir,s)
    write.table(enrichRes.ls[[s]], fp, sep = ', ', row.names = FALSE)
    message(paste("[INFO] Save to", fp))
})

## Plot top-n pvalue term for each cluster
plotGOTopn = function(go.obj, topn=20, title=NULL, fp=NULL){
    if(is.null(title)) title = "temp"
    pathway = data.frame(go.obj, stringsAsFactors = F)
    pathway = pathway[!is.na(pathway$ID),]
    if(nrow(pathway) > 20) pathway = pathway[1:topn,]
    pathway = pathway[order(pathway$pvalue, decreasing = T),]
    tt = factor(pathway$Description, levels = unique(pathway$Description))
    pp = ggplot(pathway, aes(-1*log10(p.adjust), tt))
    pbubble = pp + geom_point(aes(size = Count, color = -1*log10(p.adjust)))+
        scale_colour_gradient(low = "blue",high = "red") +
        theme_bw() +scale_size(range = c(2, 10)) +
        labs(y = "", title = title)+
        scale_y_discrete(position = 'right')+
        xlim(min(-1*log10(pathway$p.adjust))*0.9, max(-1*log10(pathway$p.adjust))*1.1)+
        theme(
            axis.text.y = element_text(size = 16, color = "black", angle = 0,vjust = 0),
            legend.position = 'none'
        )
    
    if(is.null(fp))
        return(pbubble)
    else if(fp=="return")
        return(pbubble)
    else if(fp=="print")
        print(pbubble)
    else if(!str_ends(fp,'.pdf'))
        return(pbubble)
    else
        ggsave(file = fp, plot = pbubble, width = 20, height = 10)
}

output.dir = "single_plot"
if(!dir.exists(output.dir)) dir.create(output.dir)
lapply(names(enrichRes.ls), function(name){
    ergo = enrichRes.ls[[name]]
    pathway = data.frame(ergo, stringsAsFactors = F)[1:20,]
    pathway = pathway[order(pathway$pvalue, decreasing = T),]
    tt = factor(pathway$Description, levels = unique(pathway$Description))
    pp = ggplot(pathway, aes(-1*log10(p.adjust), tt))
    pbubble = pp + geom_point(aes(size = Count, color = -1*log10(p.adjust)))+
        scale_colour_gradient(low = "blue",high = "red") +
        theme_bw() +scale_size(range = c(2, 10)) +
        labs(title = name)+
        ylab("") + 
        xlim(min(-1*log10(pathway$p.adjust))*0.9, max(-1*log10(pathway$p.adjust))*1.1)+
        theme(axis.text.y = element_text(size = 20, family = "Helvetica", color = "black", angle = 0),
              axis.text.x = element_text(size = 15, family = "Helvetica", color = "black", angle = 0))+
        theme(legend.title = element_text(color = "black", size = 20, family = "Helvetica"))+
        theme(legend.text = element_text(color="azure4", size = 15, family = "Helvetica"))+
        theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
        theme(legend.key.size = unit(1.1,'cm'))
    fp = paste0(plot_dir, '/', name, '.pdf')
    message(paste('Save to', fp, '...'))
    ggsave(file = fp, plot = pbubble, width = 20, height = 10)
})