## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"

library(reshape2)
library(cowplot)
library(ggrepel)
library(stringr)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
str2vec = function(str, sep=','){
    .split = function(s) str_split(s, sep)[[1]]
    if(length(str)>1) lapply(str, .split)
    else .split(str)
}

eve.df = readr::read_tsv(paste0(PROJECT_DIR, "/Reference/Mmus38.txt"))
eve.df = eve.df[!is.na(eve.df$chr),]
eve.df$chr = paste0("chr", eve.df$chr)
eve.df_element.type = sapply(eve.df$Repbase, function(term){
    lst = str_split(term, ' ')[[1]]
    vs = sapply(lst, function(s)str_match(s, "(\\S+)\\|(\\S+)/(\\S+)")[,3])
    vs = na.omit(vs)
    paste(vs, collapse = ',')})
eve.gr = GRanges(
    seqnames = eve.df$chr, 
    ranges = IRanges(eve.df$start, eve.df$end),
    strand = eve.df$strand)
eve.df$element.type = eve.df_element.type
eve.gr@metadata$element.type = eve.df_element.type


## Stacked barplot & distance to TSS based on peaks called by macs2
{
    narrowPeaks = readRDS(paste0(PROJECT_DIR, "/05 peak calling/narrowPeaks.rds"))
    peak.grs = lapply(narrowPeaks, function(df){
        GRanges(
            seqnames = df$chr,
            ranges = IRanges(df$start, df$end),
            strand = rep('+', nrow(df))
        )
    })
    
    peak2eve.ovs = lapply(peak.grs, findOverlaps, subject = eve.gr, maxgap = 50)
    peak2eve.sid = lapply(peak2eve.ovs, subjectHits)
    peak2eve.eve = lapply(peak2eve.sid, function(ids)eve.df[ids,])
    
    df = Reduce(rbind, peak2eve.eve)
    df = df[match(unique(df$ID), df$ID),]
    
    element.types.only = c("LTR", "LINE", "DNA", "SINE", "Unknown", "others")
    element.types = sapply(peak2eve.eve, function(df){
        x = table(Reduce(c, str2vec(df$element.type)))
        for(et in element.types.only) if(is.na(x[et])) x[et] = 0
        return(x)
    })
    element.types = element.types[element.types.only,]
    element.types = element.types[rowSums(element.types) != 0,]
    element.types = as.data.frame(t(element.types))
    
    # element.types was used for stacked barplot
    local({
        width.df = element.types
        width.df = width.df[order(rowSums(width.df)),]
        ggdf = data.frame(
            XN = rep(nrow(width.df):1, each = ncol(width.df)),
            XS = rep(rev(rownames(width.df)), each = ncol(width.df)),
            L = factor(rep(rev(colnames(width.df)), nrow(width.df)), levels = rev(colnames(width.df))),
            Y = rev(c(t(width.df))), stringsAsFactors = F
        )
        ggdf$Y = ifelse(ggdf$Y==0,NA,ggdf$Y)
        x_labels = unique(ggdf$XS)[order(unique(ggdf$XN))]
        ggplot(ggdf, aes(x = XN, y = Y, fill = L)) + 
            geom_bar(stat = 'identity', position = 'stack') +
            coord_flip() + 
            labs(x = '', y = '', title = '') +
            guides(fill = guide_legend(reverse = TRUE)) +
            theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                legend.position = 'right',
                legend.title = element_blank()
            )
        ggsave("stacked_barplot.pdf", width = 10, height = 8)
    })
}


## distance to TSS
{
    df = Reduce(rbind, peak2eve.eve)
    df = df[match(unique(df$ID), df$ID),]
    gr = GRanges(df$chr, IRanges(df$start, df$end), df$strand)
    df.anno = annotatePeak(
        peak = gr,
        TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
        tssRegion = c(-3000, 3000), 
        addFlankGeneInfo = TRUE, 
        flankDistance = 5000,
        annoDb = "org.Mm.eg.db",
        verbose = FALSE
    )
    
    X.df = local({
        X = df.anno@anno@elementMetadata@listData$distanceToTSS
        X = X[-2.5e05<X&X<2.5e05]
        X = cut(X, breaks = 40)
        X.breaks = str_match(levels(X), "\\((.*),(.*)]")[,2:3]
        X.breaks = sort(unique(as.numeric(X.breaks)))
        X.means = X.breaks[seq(1,length(X.breaks)-1)]+
            X.breaks[seq(2,length(X.breaks))]
        X.means = sapply(X.means, '/', 2)
        X.counts = `names<-`(as.integer(table(X)), NULL)
        X.ratio = X.counts / sum(X.counts)
        X.df = data.frame(breaks=X.means, counts=X.counts,ratio=X.ratio)
    }) # plot
    
    ggplot(X.df, aes(breaks, ratio)) + 
        geom_line(color='#5fb1b5', size=1.5) +
        labs(x="Absolute distance to TSS", y="EVEs proportion") + 
        ggtitle("Absolute distance to TSS of filtered EVEs") +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            axis.line = element_line(size = 0.5, colour = 'black'),
            plot.title = element_text(hjust = 0.5)
        )
    .opdir = "distance_to_TSS"
    ggsave(paste0(.opdir, ".pdf"), width = 10, height = 6)
    saveRDS(X.df, file = paste0(.opdir, ".rds"))
    write.csv(X.df, paste0(.opdir, ".csv"))
}


## Number of EVEs types in each cluster
{
    eve.types = lapply(peak2eve.eve, function(evedf){
        sapply(evedf$Repbase, function(string){
            lst = str_split(string, ' ')[[1]]
            vs = sapply(lst, function(s){
                str_match(s, "(\\S+?)[_-].*\\|(\\S+)/(\\S+)")[,1+1]
            })
            vs = na.omit(vs)
            paste(vs, collapse = ',')
        })
    })
    x = lapply(eve.types, function(et) Reduce(c, lapply(et, str2vec)))
    eve.types.keys = sort(unique(Reduce(c, x)))
    x = sapply(x, function(xc){
        y = table(xc)
        for(k in eve.types.keys)
            if(is.na(y[k]))
                y[k] = 0
            return(y[eve.types.keys])
    })
    
    type.df = as.data.frame(x)
    type.dist = dist(type.df, method = "euclidean")
    type.hc = hclust(type.dist, method = "complete")
    plot(type.hc, hang=-1)
    
    cluster.df = as.data.frame(t(x))
    cluster.dist = dist(cluster.df, method = "euclidean")
    cluster.hc = hclust(cluster.dist, method = "complete")
    plot(cluster.hc, hang=-1)
    
    x.long = `colnames<-`(melt(x), c("type","cluster","counts"))
    x.long$enrichment = scale(x.long$counts)
    x.long.type.levels = type.hc$labels[type.hc$order]
    x.long$type = factor(x.long$type, levels = x.long.type.levels)
    x.long.cluster.levels = cluster.hc$labels[cluster.hc$order]
    x.long$cluster = factor(x.long$cluster, levels = x.long.cluster.levels)
    
    cinfo <- readRDS("E:/project/MMOB/mmob_cinfo.rds")
    ggplot(x.long, aes(cluster, type))+
        geom_tile(aes(fill = enrichment), color = "white")+
        labs(x = "", y = "")+
        ggtitle("EVEs types in each cluter")+
        scale_fill_gradient(low = "white", high = "#3500a1")+
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank()
        )
    .od = "EVEs_in_cluster."
    ggsave(paste0(.od, "pdf"), width = 12, height = 10)
}


## Functional analysis of cluster-specific EVEs
{
    grs = lapply(peak2eve.eve, function(df)GRanges(df$chr, IRanges(df$start, df$end), df$strand))
    peak2eve.anno = local({
        ids = names(grs)
        names(ids) = ids
        lapply(ids, function(id){
            gr = grs[[id]]
            message(paste("[INFO] Run cluster", id))
            res = try({
                annotatePeak(
                    peak = gr,
                    TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                    tssRegion = c(-3000, 3000), 
                    addFlankGeneInfo = TRUE, 
                    flankDistance = 5000,
                    annoDb = "org.Mm.eg.db",
                    verbose = FALSE
                )
            })
            if("try-error" %in% class(res)){
                message(paste("[INFO] Error on cluster", id))
                return(NA)
            }
            else return(res)
        })
    })
    saveRDS(peak2eve.anno, file = "EVEAnnoByCluster.rds")
    peak2eve.anno = peak2eve.anno[!is.na(peak2eve.anno)]
    peak2eve.anno.genic = lapply(peak2eve.anno, function(x){
        ids = x@detailGenomicAnnotation$genic
        gs = x@anno@elementMetadata@listData$SYMBOL[ids]
        gs[gs != 'a']
    })
    peak2eve.anno.genic = peak2eve.anno.genic[sapply(peak2eve.anno.genic, length)>0]
    
    library(clusterProfiler)
    ezid.ls = lapply(peak2eve.anno.genic, function(gs){
        df = bitr(gs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        as.vector(df$ENTREZID, "character")
    })
    
    evego = lapply(ezid.ls, function(ids) enrichGO(ids, OrgDb = org.Mm.eg.db, ont = 'BP', pAdjustMethod = 'BH', readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1))
    evego = evego[!sapply(evego, is.null)]
    saveRDS(evego, file = ".tmp_evego.rds")
    if(!dir.exists(enrichRes.dir <- "evego")) 
        dir.create(enrichRes.dir)
    .=lapply(names(evego), function(s){
        fp = sprintf("%s/GO-enrich_cluster%s.csv",enrichRes.dir,s)
        write.table(evego[[s]], fp, sep = ', ', row.names = FALSE)
        message(paste("[INFO] Save to", fp))
    })
    
    #evego = readRDS(".tmp_evego.rds")
    if(!dir.exists(.opdir <- "evego_single")) dir.create(.opdir)
    cluster.sel = as.character(c(2,4,8,10,11,13,17,18))
    names(cluster.sel) = cluster.sel
    ps = lapply(cluster.sel, function(nm){
        title = paste("Cluster", nm, "-", cinfo$cell.type[cinfo$ordered.num==as.numeric(nm)])
        fp = sprintf("%s/%s.pdf", .opdir, title)
        plotGOTopn(evego[[nm]], title = title, fp = fp)
        #plotGOTopn(evego[[nm]], title = nm, fp = "return")
    })
    
    evego.dfs = lapply(evego, data.frame, stringsAsFactors = FALSE)
    evego.dfs = lapply(evego.dfs, function(df) df[!is.na(df$ID),])
    evego.df = merge2godf(evego.dfs, merge.field = "Cluster")
    inkw = c('neural','neuron','axon','spinal','brain','eye','synaptic','sensory','visual','synaps','immune','interleukin','pattern','signal','pathway','development','morphogenesis','differentiation')
    x = make_ggdf(evego.df, include = inkw, include_mode = 'reg', cinfo = cinfo)
    x = evego.df[query_idx(inkw, 'reg', godf = evego.df),]
    
    evego.dfs.topn = lapply(evego.dfs, function(x)x[1:8,])
    evego.dfs.topn.desc = lapply(evego.dfs.topn, function(x)x$Description)
    evego.dfs.topn.desc.keys = unique(Reduce(c, evego.dfs.topn.desc))
    names(evego.dfs.topn.desc.keys) = evego.dfs.topn.desc.keys
    evego.df = sapply(evego.dfs.topn, function(df){
        sapply(evego.dfs.topn.desc.keys, function(k){
            v = df$p.adjust[df$Description == k]
            if(length(v)==0) v = 0
            return(v)
        })
    })
    evego.df2 = melt(evego.df)
    colnames(evego.df2) = c("desc","cluster","value")
    .hc = hclust(dist(evego.df, method = "euclidean"), method = "complete")
    evego.df2$desc = factor(evego.df2$desc, levels = .hc$labels[.hc$order])
    .df = data.frame(t(evego.df), stringsAsFactors = F)
    .hc = hclust(dist(.df, method = "euclidean"), method = "complete")
    evego.df2$cluster = factor(evego.df2$cluster, levels = .hc$labels[.hc$order])
    p = {
        ggplot(evego.df2, aes(cluster, desc))+geom_tile(aes(fill = value), color = "white")+
            labs(x = "Cluster", y = "")+ggtitle("Top8 functional enrichment of each cluster")+
            scale_y_discrete(position = "right")+
            scale_fill_gradient(low = "white", high = "#3500a1")+
            theme(
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA),
                plot.title = element_text(size = 20, hjust = 0.5),
                axis.ticks = element_blank(),
                axis.text.y = element_text(size = 16),
                axis.text.x = element_text(size = 16),
                legend.position = "left",
            )
    }
    ggsave(".tmp_evego_topnplot.pdf",p,width = 14, height = 16)
    saveRDS(p, ".tmp_evego_topnplot.rds")
}


## Venn plot
{
    # Use `peak2eve.sid`
    ids = unique(cinfo$cell.type)
    eve.in.ct = lapply(`names<-`(ids,ids), function(celltype){
        cl.sel = cinfo$ordered.num[cinfo$cell.type == celltype]
        unique(Reduce(c, peak2eve.sid[cl.sel]))
    })
    #library(VennDiagram)
    #venn.diagram(eve.in.ct, filename = NULL)
    .output_dir = "eve.in.ct"
    if(!dir.exists(.output_dir)) dir.create(.output_dir)
    .=lapply(names(eve.in.ct), function(s){
        fp = sprintf("%s/%s.csv", .output_dir, s)
        df = data.frame(paste(eve.in.ct[[s]]))
        write.table(df,fp,col.names=F,row.names=F,quote=F,sep='\n')
    })
    
    .df = data.frame(
        ct = names(eve.in.ct),
        id = sapply(eve.in.ct, function(ids)paste(ids, collapse = ','))
    )
    write.table(.df, "eve.in.ct.csv",quote = F,col.names = F,row.names = F)
}

save.image("EVE.RDATA")