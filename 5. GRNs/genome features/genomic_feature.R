library(ChIPseeker)

cinfo = readRDS(paste0(PROJECT_DIR, "/04 Clustering/02 annotate celltypes/cinfo.rds"))
cinfo = `rownames<-`(cinfo[order(cinfo$cell.type),], NULL)
cinfo = `rownames<-`(cinfo[c(6:14,21,20,4,5,1,16:19,15,2,3),], NULL)

peakAnno = readRDS(paste0(PROJECT_DIR, "/05 peak calling/peakAnno.ls.rds"))
peakAnno2 = `names<-`(peakAnno, paste0('C', names(peakAnno)))
peakAnno2 = peakAnno2[paste0('C', cinfo$ordered.num)]
names(peakAnno2) = paste0(names(peakAnno2), '(', cinfo$cell.type, ')')
plotAnnoBar(peakAnno2)
ggsave("GenomicFeature.pdf", width = 12, height = 6)
