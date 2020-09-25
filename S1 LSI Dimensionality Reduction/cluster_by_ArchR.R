suppressMessages(library(ArchR))
addArchRGenome("mm10")
addArchRThreads(threads = 8)

mob = ArchRProject("/home/xiangrong1/hw5_COHORT/project/mmob/ArchRAnalysis/mob.arrow", "output", FALSE)

## intersect cells
snap_barcodes = readRDS("/home/xiangrong1/hw5_COHORT/project/mmob/link/xsp__barcode.rds")
arch_barcodes = sapply(rownames(mob), function(s){strsplit(s,'#')[[1]][2]})
intersect_bcs = intersect(snap_barcodes, arch_barcodes) # 0.9974407
mob = mob[paste0('mob#', intersect_bcs),]

## independent clustering
mob = addIterativeLSI(
    ArchRProj = mob,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list(
        resolution = c(1), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

mob = addUMAP(
    ArchRProj = mob, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

mob <- addTSNE(
    ArchRProj = mob, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30,
    force = TRUE
)

mob <- addClusters(
    input = mob,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "ArchrClusters2",
    resolution = 2,
    force = TRUE
)

table(mob$ArchrClusters)

# add snap cluster to mob
snap_clusters = readRDS("/home/xiangrong1/hw5_COHORT/project/mmob/link/xsp__cluster.rds")
names(snap_clusters) = snap_barcodes
snap_clusters = snap_clusters[intersect_bcs]
cinfo = readRDS("/home/xiangrong1/hw5_COHORT/project/mmob/link/cinfo.rds")
cinfo$cellTypes = c('OEC','End1','EX2','AST','OLI1','MG','OSN','OPC','Purk','IN5','End1','IN7','IN8','Gran','IN4','EX1','IN3','IN1','IN2','IN6','OLI2')
snap_clusters = cinfo$cellTypes[match(snap_clusters, cinfo$ordered.num)]
mob$SnapClusters = snap_clusters

#cM <- confusionMatrix(mob$ArchrClusters, mob$SnapClusters)
cM <- confusionMatrix(mob$SnapClusters, mob$ArchrClusters)
saveRDS(cM, file = "confusionMatrix.rds")
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black",
    filename = "confusionMatrix-ArchR-SnapATAC-clusters2.pdf"
)

## reorder
cnames = c('IN1','IN2','IN8','IN3','IN4','IN5','IN6','IN7','EX1','EX2','OSN','Purk','OEC','Gran','AST','OLI1','OLI2','OPC','MG','End1')
rnames = c('C17','C16','C15','C9','C14','C12','C21','C19','C20','C10','C13','C3','C4','C11','C18','C6','C7','C1','C2','C8','C5')
#cM2 = as.matrix(cM)[rnames, cnames]
cM2 = as.matrix(cM)[cnames, rnames]
pheatmap::pheatmap(
    mat = cM2, 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black",
    cluster_rows = F,
    cluster_cols = F,
    filename = "confusionMatrix-ArchR-SnapATAC-clusters2.pdf"
)

p1 = plotEmbedding(ArchRProj = mob, colorBy = "cellColData", name = "ArchrClusters", embedding = "UMAP")
p2 = plotEmbedding(ArchRProj = mob, colorBy = "cellColData", name = "SnapClusters", embedding = "UMAP")
plotPDF(p1, p2, name = "Plot-UMAP-ArchR-SnapATAC.pdf", ArchRProj = mob, addDOC = FALSE, width = 6, height = 6)


## plot marker heatmap

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
    #markers = filterMarkers(markers, ref)
    return(markers)
}
#markers.df = load_markers("/home/xiangrong1/hw5_COHORT/project/mmob/link/markers.txt")
markers.df = load_markers("/home/xiangrong1/hw5_COHORT/project/mmob/ArchRAnalysis/ClusteringComparation/markers")
plotPDF(
    plotMarkerHeatmap(
        seMarker = markersGS, 
        cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
        labelMarkers = markers.df$symbol,
        transpose = FALSE
    ),
    name = paste0("MarkerHeatmap-ArchR2"), 
    width = 6, height = 10, ArchRProj = mob, addDOC = FALSE
)


mob$ArchrClusters = as.character(mob$ArchrClusters)
mob$SnapClusters = as.character(mob$SnapClusters)

local({
    for (ctype in unique(markers.df$ctype)) {
        markers = markers.df$symbol[markers.df$ctype == ctype]
        markersGS = getMarkerFeatures(
            ArchRProj = mob,
            useMatrix = "GeneScoreMatrix",
            groupBy = "ArchrClusters",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "wilcoxon"
        )
        plotPDF(
            plotMarkerHeatmap(
                seMarker = markersGS, 
                cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
                labelMarkers = markers,
                transpose = TRUE
            ),
            name = paste0("MarkerHeatmap-", ctype), 
            width = 8, height = 6, ArchRProj = mob, addDOC = FALSE
        )
    }
}) # plot markers heatmap for each cell type

saveArchRProject(ArchRProj = mob, outputDirectory = "ClusteringComparation", load = FALSE)