## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"

library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

peak = readRDS(paste0(PROJECT_DIR, "/xsp__peak.rds"))
idy.ls = readRDS(paste0(PROJECT_DIR, "/idy.ls.rds"))
idp.ls = lapply(idy.ls, function(vec){peak[vec]})

## Run ChIPseeker::annotatePeak for each cluster
peakAnno.ls = lapply(idp.ls, function(gr){
    annotatePeak(gr, TxDb=txdb, verbose=FALSE, tssRegion=c(-3000, 3000), 
                 addFlankGeneInfo=TRUE, flankDistance=5000, 
                 annoDb="org.Mm.eg.db")
})
saveRDS(peakAnno.ls, file = "peakAnno.ls.rds")