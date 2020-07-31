library(SnapATAC)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
xsp = readRDS(paste0(PROJECT_DIR, "/xsp.rds"))

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
saveRDS(DARS, "DARs.rds")

idy.ls = lapply(DARS, function(d){
    d$FDR = p.adjust(d$PValue, method="BH")
    idy = which(d$FDR < 5e-2 & d$logFC > 0)
    if(length(idy) > 1000){
        PValues = d$PValue
        PValues[d$logFC < 0] = 1
        idy = order(PValues, decreasing = FALSE)[1:1000]
    }
    idy
})

saveRDS(idy.ls, file = paste0(PROJECT_DIR, "/idy.ls.rds"))