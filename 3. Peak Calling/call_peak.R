library(SnapATAC)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
xsp = readRDS(paste0(PROJECT_DIR, "/xsp.rds"))

output_dir_macs2 = sprintf("%s/output_macs2", getwd())
if(!dir.exists(output_dir_macs2)) dir.create(output_dir_macs2)

setwd(output_dir_macs2)
path2snaptools = "/path/snaptools"
path2macs = "/path/macs2"
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
    
    narrowPeaks = lapply(fps, function(fp){
        df = read.table(fp, sep = '\t', header = F, stringsAsFactors = F)
        colnames(df) = c('chr','start','end','name','intm10log10qvalue','strand','fold.change','mlog10pvalue','mlog10qvalue','relative.summit.position.to.peak.start')
        return(df)
    })
    
    ## Set pattern according to your file name
    ## Use () to catch the keyword presenting cluster
    filenamePattern = 'C(\\d+)\\.narrowPeak'
    
    names(narrowPeaks) = stringr::str_match(basename(fps), filenamePattern)[,2]
    narrowPeaks = narrowPeaks[order(as.integer(names(narrowPeaks)))]
    
    saveRDS(narrowPeaks, file = "narrowPeaks.rds")
})

## Add pmat to xsp
xsp = createPmat(xsp, peak.gr, 20, TRUE, num.cores = 4)

saveRDS(xsp@peak, file=paste0(PROJECT_DIR, "/xsp__peak.rds"))
saveRDS(xsp, file=paste0(PROJECT_DIR, "/xsp.rds"))