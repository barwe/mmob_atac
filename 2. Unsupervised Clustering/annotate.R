library(SnapATAC)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
xsp = readRDS(paste0(PROJECT_DIR, "/xsp.rds"))


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


gene.bed = readRDS(paste0(PROJECT_DIR, "/mm10.geneBed.rds"))
markers = load_markers('markers.txt', gene.bed$symbol)
gr = GRanges(gene.bed$chr, IRanges(gene.bed$start, gene.bed$end))
xsp = create_gmat(xsp, markers, gr)
gmat = save_gmat(xsp, output.filename = paste0(PROJECT_DIR, "/xsp__gmat.rds"))
save_single_feature_tsne(xsp, markers, gmat)

saveRDS(xsp, file=paste0(PROJECT_DIR, "/xsp.rds"))