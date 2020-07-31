## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
PATH_TO_HOMER = "/path/homer2"

xsp = readRDS("/home/xiangrong1/mmob/xsp.rds")
cbp = readRDS("/home/xiangrong1/mmob/cell_peak_dgCMatrix.rds")
xsp@pmat = cbp; range(cbp)
#idy.ls = readRDS("/home/xiangrong1/hw5/project/mmob/idy.ls.0519.rds")

motifs = list()
path.to.homer = system("which homer", intern = T)
for (cidx in 1:21){
    motifs[cidx] = runHomer(
        xsp[, idy.ls[[cidx]], 'pmat'],
        mat = 'pmat',
        path.to.homer = "/zfssz2/ST_MCHRI/COHORT/wangshiyou/software/Homer/bin/findMotifsGenome.pl",
        result.dir = paste("motifs/", cidx, sep = ''),
        num.cores = 5,
        genome = 'mm10',
        motif.length = 10,
        scan.size = 300,
        optimize.count = 2,
        background = 'automatic',
        local.background = FALSE,
        only.known = FALSE,
        only.denovo = FALSE,
        fdr.num = 5,
        cache = 100,
        overwrite = TRUE,
        keep.minimal = FALSE
    )
}

motif.ls = lapply(1:21, function(i){
    fp = sprintf("motifs/%s/knownResults.txt", i)
    df = read.table(fp, sep = '\t', header = T, comment.char = '')
    colnames(df) = c(colnames(df)[1:5], 'target.count', 'target.freq', 'background.count', 'background.freq')
    col_names = colnames(df)
    df$Cluster = rep(i, nrow(df))
    df = df[,c('Cluster', col_names)]
    return(df)
})

motif = Reduce(rbind, motif.ls)

saveRDS(motif.ls, file = "motif.ls.rds")
saveRDS(motif, file = "motif.rds")
