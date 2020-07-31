library(SnapATAC)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"

make.speaks = function(symbols, spmat = 21){
    # Local function `eval_makrer` need `idy.ls` data
    eval_marker = function(marker, idy = idy.ls, clusters = NULL){
        if(is.null(clusters)) clusters = names(idy)
        marker_pos = c(marker$start, marker$end)
        scores = `names<-`(sapply(clusters, function(cx){
            da.peaks = peak.bed[idy[[cx]],]
            validate.by.chr = which(as.character(da.peaks$chr) == as.character(marker$chr))
            peak_pos =da.peaks[validate.by.chr, c('start', 'end')]
            isolp = isOverlapping(peak_pos, marker_pos, return.logical = TRUE)
            sc = sum(isolp)
        }), clusters)
        return(scores)
    }
    isOverlapping = function(target, ref, orate = 0, return.logical = TRUE){
        #`is.ratio: return logical vector or numeric overlapping ratio
        N = nrow(target)
        target.len = target[[2]] - target[[1]]
        ref.len = ref[2] - ref[1]
        x.len = target.len + ref.len
        x.min = sapply(target[[1]], function(x) min(x, ref[1]))
        x.max = sapply(target[[2]], function(x) max(x, ref[2]))
        x.range = x.max - x.min
        if(return.logical)
            isoverlapping = x.range < x.len
        else{
            isoverlapping = (x.len - x.range) / target.len
            isoverlapping = ifelse(x.len - x.range == ref.len, 1, isoverlapping)
            isoverlapping = ifelse(isoverlapping <= orate, 0, isoverlapping)
        }
        return(isoverlapping)
    }
    #----------------#
    if(is.numeric(spmat)){
        row_names = 1:spmat
        spmat = data.frame(row.names = row_names)
    }
    else {
        row_names = rownames(spmat)
        symbols.ignored = which(symbols %in% colnames(spmat))
        if(length(symbols.ignored) > 0){
            .tmp = paste(symbols[symbols.ignored], collapse = ', ')
            message(paste('Ignore markers:', .tmp))  
            symbols = symbols[-symbols.ignored]
        }
    }
    msg = list(zero = c(NULL), moreThanOne = c(NULL))
    pb = progress::progress_bar$new(total = length(symbols))
    res = `names<-`(lapply(symbols, function(symbol){
        pb$tick()
        marker = gene.bed[which(gene.bed$symbol == symbol),]
        if(nrow(marker) == 0) {
            msg$zero <<- c(msg$zero, symbol)
            return(rep(0, length(row_names)))
        } else if(nrow(marker) == 1) {
            res = eval_marker(marker, clusters = row_names)
            return(res)
        } else {
            msg$moreThanOne <<- c(msg$moreThanOne, symbol)
            res = eval_marker(marker[1,], clusters = row_names)
            return(res)
        }
    }), symbols)
    tmp.fn = function(x){
        tmp.str = paste0(x, ': ', paste(msg[[x]], collapse = ', '))
        message(tmp.str)
    }
    lapply(names(msg), tmp.fn)
    spmat = as.data.frame(res, row.names = row_names)
    return(spmat[,symbols,drop=FALSE])
}
plotdot = function(speak, ctypes = NULL, ylabs = NULL, symbols = NULL, clusters = NULL, scale.pt = 1, title = '', hc.order = NULL ){
    ##################################################################
    ##`speak: a density dataframe with clusters row and markers column
    ##`symbols: a markers vector, default `colnames(speak)`
    ##`clusters: a clusters vector, default `rownames(speak)`
    ##`scale.pt: scale the point size, default 1
    ##`title: title of figure
    ##`ylabs: labels shown as y axis, default `colnames(speak)`
    ##`hc.order: cluster index corresponding with hc
    ##################################################################
    library(ggplot2)
    if(is.null(symbols)) symbols = colnames(speak)
    if(is.null(clusters)) clusters = rownames(speak)
    if(is.null(ylabs)) ylabs = symbols
    if(is.null(hc.order)) hc.order = clusters
    xlabels = hc.order[which(hc.order %in% clusters)]
    if(is.null(ctypes)) xlabs = xlabels
    else{
        cts = sapply(xlabels, function(x){
            v = ctypes[[paste0('C', x)]]
            if(is.null(v)) return('')
            else return(v)
        })
        xlabs = paste(xlabels, cts, sep = '\n')
    }
    numClusters = nrow(speak)
    numMarkers = ncol(speak)
    #-----------------------
    speak = speak[clusters,,drop = FALSE][,rev(symbols),drop = FALSE]
    speak = speak[xlabels,,drop = FALSE]
    pt.size = Reduce(c, speak)
    pt.alpha = unique(pt.size)
    pt.alpha = match(pt.size, sort(pt.alpha)) / max(pt.alpha)
    df = data.frame(
        Cluster = rep(seq(numClusters), numMarkers),
        Marker = rep(seq(numMarkers), each = numClusters),
        PointSize = Reduce(c, speak)^scale.pt,
        PointColor = ifelse(pt.size == 0, 'F', 'T'),
        PointLabel = ifelse(pt.size == 0, '', pt.size),
        PointAlpha = pt.alpha + 0.5
    )
    p = ggplot(df, aes(Cluster, Marker))+ 
        geom_point(aes(size = PointSize, color = PointColor, alpha = PointAlpha))+ 
        labs(x = "", y = "", title = title)+ 
        scale_x_continuous(breaks = seq(numClusters), labels = xlabs)+ 
        scale_y_continuous(breaks = seq(numMarkers), labels = rev(ylabs))+ 
        scale_color_manual(values = c('white', 'red'), limits = c('F', 'T'))+ 
        scale_size(range = c(1, 4)) +
        theme(
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "top",
            axis.text.x = element_text(size = 6)
        )
}

idy.ls = readRDS(paste0(PROJECT_DIR, "/idy.ls.rds"))
gene.bed = readRDS(paste0(PROJECT_DIR, "/Reference/mm10.geneBed.rds"))
peak.bed = local({
    peak.gr = readRDS("PeakCalling/xsp@peak.rds")
    chr = as.character(peak.gr@seqnames)
    start = peak.gr@ranges@start
    end = peak.gr@ranges@start + peak.gr@ranges@width - 1
    data.frame(chr=chr,start=start,end=end,stringsAsFactors = F)
})
markers = load_markers('ClusterAnno/markers.txt', gene.bed$symbol)

xsp_cluster = readRDS("xsp_cluster.rds")
speaks = make.speaks(markers$symbol, spmat = length(levels(xsp_cluster)))
cinfo = readRDS("ClusterAnno/cluster_anno.rds")
speaks = speaks[cinfo$ordered.num, markers$symbol]

p = plotdot(speaks, paste0(markers$symbol, ' (', markers$ctype, ')'))
ggsave("dotplot.pdf",width=6,height=10,limitsize=FALSE)