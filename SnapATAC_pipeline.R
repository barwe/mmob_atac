library(GenomicRanges)
library(SnapATAC)
library(leiden)
library(IRanges)
library(GenomicRanges)
library(stringr)
library(ggplot2)

##========================================================##
##                 Step01: Select barcode                 ##
##=======================================================

#@ Create SnapATAC object from prepared .snap database
xsp = createSnap(
    file = 'mouse_olfactory_bulb_mm10_2012.snap',
    sample = 'mouse_olfactory_bulb_mm10_2012'
)

#@ Plot barocde quality control distribution
plotBarcode(xsp, "plot01-barcode_info.pdf",
            pdf.width = 7, pdf.height = 7,
            col = 'grey', border = 'grey',
            breaks = 50)

#@ Filter cells based on fragment.num and UMI
## fragment.num - total number of fragments per barcode
## UMI - unique molecular identifier
xsp = filterCells(
    obj = xsp,
    subset.names = c('fragment.num', 'UMI'),
    low.thresholds = c(1000,10),
    high.thresholds = c(Inf, Inf)
)

##----------------------------------------------------------



##========================================================##
##                  Step02: Create bmat                   ##
##=======================================================

showBinSizes("mouse_olfactory_bulb_mm10_2012.snap")
xsp = addBmatToSnap(xsp, bin.size=5000, num.cores=1)

##----------------------------------------------------------



##========================================================##
##                  Step03: Select bins                   ##
##=======================================================

#@ Binarize the bmat
xsp = makeBinary(xsp, mat="bmat")

#@ Remove bins in ENCODE blacklist
blacklist = readRDS('mm10.blacklist.rds')
blacklist.ir = IRanges(blacklist$start, blacklist$end)
blacklist.gr = GRanges(blacklist$chr, blacklist.ir)
idx_rm = queryHits(findOverlaps(xsp@feature, blacklist.gr))
if(length(idx_rm) > 0) xsp = xsp[, -idx_rm, mat="bmat"]
rm(blacklist, blacklist.ir, blacklist.gr, idx_rm)

#@ Remove bins located at unwanted chromosomes
grep_idx = grep("random|chrM|chrUn", seqlevels(xsp@feature))
chr_rm = seqlevels(xsp@feature)[grep_idx]
idx_rm = grep(paste(chr_rm, collapse="|"), xsp@feature)
if(length(idx_rm) > 0) xsp = xsp[, -idx_rm, mat="bmat"]
rm(grep_idx, chr_rm, idx_rm)

#@ Cut off bins
bin.cov = log10(1 + Matrix::colSums(xsp@bmat))
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idx_sel = which(bin.cov <= bin.cutoff & bin.cov > 0)
if(length(idx_sel) > 0) xsp = xsp[, idx_sel, mat="bmat"]
rm(bin.cov, bin.cutoff, idx_sel)

##----------------------------------------------------------



##========================================================##
##            Step04: Select principal component          ##
##=======================================================

xsp = runDiffusionMaps(xsp, "bmat", num.eigs=50)
#saveRDS(xsp@jmat, file="xsp_jmat.rds")
#saveRDS(xsp@smat, file="xsp_smat.rds")
#saveRDS(xsp@regModel, file="xsp_regModel.rds")

plotDimReductPW(obj=xsp, eigs.dims=1:50, down.sample=5000,
                point.size=0.3, point.color="grey", 
                point.shape=19, point.alpha=0.6,
                pdf.height=7, pdf.width=7, 
                pdf.file.name="plot04-dim_reduct.pdf")

##----------------------------------------------------------



##========================================================##
##                  Step05: Run cluster                   ##
##=======================================================

xsp = runKNN(xsp, eigs.dims = 1:20, k = 15)

xsp = runCluster(
    obj = xsp,
    tmp.folder = tempdir(),
    louvain.lib = 'leiden',
    resolution = 0.5,
    seed.use=10
)

print(length(levels(xsp@cluster)))

##----------------------------------------------------------



##========================================================##
##                    Step06: Visualize                   ##
##=======================================================

.run_viz = function(obj, type = "tsne"){
    if(type == "tsne") type = "Rtsne"
    runViz(
        obj = obj, 
        tmp.folder = tempdir(),
        dims = 2,
        eigs.dims = 1:20, 
        method = type,
        seed.use = 10
    )
}

.plot_viz = function(obj, type = "tsne"){
    title = "Mouse Olfactory Bulb scATAC-seq"
    output.filename = sprintf("plot06-%s.pdf", type)
    plotViz(
        obj = xsp,
        method = type, 
        main = title,
        point.color = xsp@cluster, 
        point.size = 0.4, 
        point.shape = 19, 
        point.alpha = 0.8, 
        text.add = TRUE,
        text.size = 1,
        text.color = "black",
        text.halo.add = TRUE,
        text.halo.color = "white",
        text.halo.width = 0.1,
        down.sample = 10000,
        legend.add = TRUE,
        legend.pos = 'topright',
        pdf.file.name = output.filename,
        pdf.width = 8,
        pdf.height = 8
    )
}

visual_method = "tsne"
xsp = .run_viz(xsp, visual_method)
.plot_viz(xsp, visual_method)

saveRDS(xsp, file="xsp.rds")

##----------------------------------------------------------



##========================================================##
##               Step07: Annotate cell type               ##
##=======================================================

## Load local markers list
load_markers = function(fp, ref){
    #`fp: marker file path
    #`ref: reference markers to filter
    loadMarkers = function(fp){
        df = read.table(fp,sep=':',stringsAsFactors=FALSE,strip.white=TRUE)
        tmp.fn = function(x)stringr::str_trim(strsplit(x, ',')[[1]])
        rt = `names<-`(lapply(df$V2, tmp.fn), df$V1)
    }
    formatMarkers = function(mks){
        tmp.fn = function(ctype) rep(ctype, length(mks[[ctype]]))
        cts = lapply(names(mks), tmp.fn)
        cts = Reduce(c, cts)
        mks = Reduce(c, mks)
        data.frame(ctype = cts, symbol = mks, stringsAsFactors = FALSE)
    }
    filterMarkers = function(df, ref){
        idx.rm = which(! df$symbol %in% ref)
        if(length(idx.rm) > 0){
            message('Remove markers not in gene.bed ...')
            print(df[idx.rm,])
        }
        return(df[which(df$symbol %in% ref),])
    }
    markers = loadMarkers(fp)
    markers = formatMarkers(markers)
    markers = filterMarkers(markers, ref)
    return(markers)
}
## Create GenomicsRanges object from local gene bed file
create_GR = function(fp){
    genes = read.table(fp)
    iranges = IRanges(genes[,2], genes[,3])
    genes.gr = GRanges(genes[,1], iranges, name=genes[,4])
}
## Create gmat based on loaded markers and GR object
create_gmat = function(obj, markers, gr, gmat=NULL){
    #@ Remove those markers whose gmat data exist if gmat is given
    if(!is.null(gmat)){
        index = markers$symbol %in% setdiff(markers$symbol, colnames(gmat))
        markers = markers[index,]        
    }
    if(nrow(markers) == 0)
        stop('No markers')
    else{
        genes.sel.gr = gr[which(gr$name %in% markers$symbol)]
        message('run SnapATAC::createGmatFromMat ...')
        obj = createGmatFromMat(
            obj = obj,
            input.mat = "bmat",
            genes = genes.sel.gr,
            do.par = TRUE
        )
        message('run SnapATAC::scaleCountMatrix ...')
        obj = scaleCountMatrix(
            obj = obj,
            cov = obj@metaData$UQ + 1,
            mat = "gmat",
            method = "RPM"
        )
        message('run SnapATAC::runMagic ...')
        obj = runMagic(
            obj = obj,
            input.mat = "gmat",
            step.size = 3
        )
    }
    return(obj)
}
## Save and update gmat matrix
save_gmat = function(obj, gmat=NULL, output.filename=NULL){
    if(is.null(output.filename))
        output.filename = paste0(getwd(), 'xsp_gmat.rds')
    if(is.null(gmat))
        gmat = obj@gmat
    else{
        gmat = cbind(gmat, obj@gmat)
        gmat = gmat[,unique(colnames(gmat))]
    }
    saveRDS(gmat, file = output.filename)
    message(paste('Save gmat to', output.filename))
    return(gmat)
}
## Plot TSNE figure for each marker based on loaded markers list
save_single_feature_tsne = function(obj, markers, gmat=NULL, output.dir=NULL){
    if(is.null(gmat))
        gmat = obj@gmat
    if(is.null(output.dir))
        output.dir = 'single_feature_tsne'
    if(!str_ends(output.dir,'/')) 
        output.dir = sprintf("%s/",output.dir)
    if(!dir.exists(output.dir)) 
        dir.create(output.dir)
    for(idx in 1:nrow(markers)){
        gene = markers$symbol[idx]
        type = markers$ctype[idx]
        fp = sprintf("%s/%s-%s.pdf", output.dir, type, gene)
        if(!file.exists(fp)){
            message(sprintf("Save to %s", fp))
            plotFeatureSingle(
                obj = obj,
                feature.value = gmat[,gene],
                method = "tsne",
                main = gene,
                point.size = 0.1,
                point.shape = 19,
                down.sample = 10000,
                quantiles = c(0, 1),
                pdf.file.name = fp
            )
        }
        else message(sprintf("Esacpe %s", fp))
    }
}

gene.bed = readRDS("mm10.geneBed.rds")
markers = load_markers('markers.txt', gene.bed$symbol)
gr = GRanges(gene.bed$chr, IRanges(gene.bed$start, gene.bed$end))
xsp = create_gmat(xsp, markers, gr)
gmat = save_gmat(xsp)
save_single_feature_tsne(xsp, markers, gmat)

##----------------------------------------------------------



##========================================================##
##          Step08: Run hierarchical clustering           ##
##=======================================================

local({
    fp = 'plot08-hclust.pdf'
    ensemble.ls = lapply(
        X = split(seq(length(xsp@cluster)), xsp@cluster),
        FUN = function(x){
            SnapATAC::colMeans(xsp[x,], mat = "bmat")
        }
    )
    hc = hclust(
        d = as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))),
        method = "ward.D2"
    )
    pdf(file = fp)
    plot(hc, hang=-1, xlab="")
    dev.off()
})

##----------------------------------------------------------



##========================================================##
##                    Step09: Call peak                   ##
##=======================================================

output_dir_macs2 = sprintf("%s/output_macs2", getwd())
if(!dir.exists(output_dir_macs2)) dir.create(output_dir_macs2)

setwd(output_dir_macs2)
path2snaptools = "/zfssz2/ST_MCHRI/COHORT/wangshiyou/software/anaconda3/bin/snaptools"
path2macs = "/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/bin/macs2"
peak.macs = runMACSForAll(
    obj = xsp,
    output.prefix = "mmob",
    path.to.snaptools = path2snaptools,
    path.to.macs = path2macs,
    gsize = "mm",
    tmp.folder = tempdir(),
    num.cores = 5,
    min.cells = 0,
    buffer.size = 500,
    macs.options = "--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
    keep.minimal = TRUE
)
setwd('..')

peak.gr = local({
    cmd = sprintf("ls %s | grep narrowPeak", output_dir_macs2)
    fps = paste0(output_dir_macs2, '/', system(cmd, intern = TRUE))
    grs = lapply(fps, function(fp){
        df = read.table(fp)
        GRanges(df[,1], IRanges(df[,2], df[,3]))
    })
    gr = reduce(Reduce(c, grs))
})

## Add pmat to xsp
xsp = createPmat(xsp, peak.gr, 20, TRUE, num.cores = 4)

## Save peaks seprately
local({
    chr = as.character(xsp@peak@seqnames)
    start = xsp@peak@ranges@start
    end = start + xsp@peak@ranges@width - 1
    df = data.frame(chr = chr, start = start, end = end)
    saveRDS(df, file = "mmob_peakBed.rds")
})

##----------------------------------------------------------



##========================================================##
##                  Step10: Identify DARs                 ##
##=======================================================

cbp = readRDS("/home/xiangrong1/mmob/cell_peak_dgCMatrix.rds")
xsp@pmat = cbp

DARS = lapply(levels(xsp@cluster), function(i){
    findDAR(
        obj = xsp,
        input.mat = "pmat",
        cluster.pos = i,
        cluster.neg = NULL,
        cluster.neg.method = "knn",
        bcv = 0.1,
        test.method = "exactTest",
        seed.use = 10
    )
})
names(DARS) = levels(xsp@cluster)

idys = (function(dars=DARS, fdr=1e-4, logfc=2, min=1000){
    idy.ls = lapply(dars, function(d){
        d$FDR = p.adjust(d$PValue, method="BH")
        idy = which(d$FDR < fdr & d$logFC > logfc)
        if(length(idy) < min){
            PValues = d$PValue
            PValues[d$logFC < 0] = 1
            idy = order(PValues, decreasing = FALSE)[1:min]
        }
        idy
    })
    print(sapply(idy.ls, length))
    return(idy.ls)
})()

saveRDS(DARS, "mmob_DARs.rds")
saveRDS(idys, "mmob_idy.ls.rds")

##----------------------------------------------------------



##========================================================##
##            Step11: Plot gene marker bubble             ##
##=======================================================
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
plotdot = function(speak, ctypes = NULL, ylabs = NULL,
   symbols = NULL, clusters = NULL, scale.pt = 1, title = '', hc.order = NULL ){
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
.timefp = function(prefix='temp', file.type='pdf', fmt=NULL){
    dir_keys = c('d','dir','directory','cat','category')
    if(is.null(fmt)) 
        fmt = 'YYYYMMDDHHMMSS'
    if(is.null(file.type)) 
        file.type = ''
    else if(nchar(file.type)==0) 
        file.type = ''
    else if(file.type%in%dir_keys)
        file.type = ''
    else
        file.type = paste0('.', file.type)
    
    time_str = as.character(Sys.time())
    ptn1 = paste(rep('(.*)',3),collapse='-')
    ptn2 = paste(rep('(.*)',3),collapse=':')
    pattern = paste(ptn1, ptn2)
    d = stringr::str_match(time_str, pattern)[1,2:7]
    names(d) = c('YYYY', 'MM', 'DD', 'HH', 'MM', 'SS')
    for (key in names(d)) {
        fmt = gsub(key, d[key], fmt, fixed = TRUE)
    }
    paste0(prefix, '_', fmt, file.type)
}

idy.ls = readRDS("DARs/differentiallyAccessiblePeakID.ls.rds")
gene.bed = readRDS("RefDatabase/mm10.geneBed.rds")
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
ggsave(.timefp("plot11-dotplot",fmt="MMDD-HHMMSS"),width=6,height=10,limitsize=FALSE)

##----------------------------------------------------------



##========================================================##
##                  Step12: Find motifs                   ##
##=======================================================

#.PATH_TO_HOMER = "/zfssz2/ST_MCHRI/COHORT/wangshiyou/software/Homer/bin/homer2"
.PATH_TO_HOMER = system("which findMotifsGenome.pl", intern = T)

#xsp@peak = readRDS("/home/xiangrong1/mmob/xsp__peak.rds")

motifs = local({
    output.dir = "motifs2/"
    motifs = list()
    if (dir.exists(output.dir)){
        system(paste("rm -rf", output.dir))
        message(paste("[INFO] Remove dir", output.dir))
    } else {
        dir.create(output.dir)
        message(paste("[INFO] Create dir", output.dir))
    }
    for (cidx in names(idys)){
        result.dir = paste0("motifs/", cidx)
        if(dir.exists(result.dir)){
            system(paste("rm -rf", result.dir))
            message(paste("[INFO] Remove dir", result.dir))
        }
        message(paste("[INFO] Run homer on cluster", cidx))
        motifs[cidx] = runHomer(
            xsp[, idys[[cidx]], 'pmat'],
            mat = 'pmat',
            path.to.homer = .PATH_TO_HOMER,
            result.dir = paste0("motifs/", cidx),
            num.cores = 5,
            genome = 'mm10',
            motif.length = 10,
            scan.size = 200,
            optimize.count = 2,
            background = 'automatic',
            local.background = FALSE,
            only.known = TRUE,
            only.denovo = FALSE,
            fdr.num = 5,
            cache = 500,
            overwrite = TRUE,
            keep.minimal = FALSE
        )
    }
    motifs
})

##----------------------------------------------------------
