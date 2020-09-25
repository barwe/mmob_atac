
## Cross-annotation based on ArchR clustering and SnapATAC clustering

IntClusters = mob$ArchrClusters
IntClusters[intersect(which(mob$ArchrClusters=="C18"),which(mob$SnapClusters=="IN2"))]="C18IN2"
IntClusters[intersect(which(mob$ArchrClusters=="C18"),which(mob$SnapClusters=="Gran"))]="C18Gran"
IntClusters[intersect(which(mob$ArchrClusters=="C17"),which(mob$SnapClusters=="IN1"))]="C17IN1"
IntClusters[intersect(which(mob$ArchrClusters=="C17"),which(mob$SnapClusters=="IN2"))]="C17IN2"
IntClusters[intersect(which(mob$ArchrClusters=="C17"),which(mob$SnapClusters=="Gran"))]="C17Gran"
IntClusters[intersect(which(mob$ArchrClusters=="C15"),which(mob$SnapClusters=="IN7"))]="C15IN7"
IntClusters[intersect(which(mob$ArchrClusters=="C15"),which(mob$SnapClusters=="IN1"))]="C15IN1"
IntClusters[intersect(which(mob$ArchrClusters=="C15"),which(mob$SnapClusters=="IN2"))]="C15IN2"
IntClusters[intersect(which(mob$ArchrClusters=="C12"),which(mob$SnapClusters=="IN5"))]="C12IN5"
IntClusters[intersect(which(mob$ArchrClusters=="C12"),which(mob$SnapClusters=="IN6"))]="C12IN6"
IntClusters[intersect(which(mob$ArchrClusters=="C21"),which(mob$SnapClusters=="IN8"))]="C21IN8"
IntClusters[intersect(which(mob$ArchrClusters=="C9"),which(mob$SnapClusters=="IN3"))]="C9IN3"
IntClusters[intersect(which(mob$ArchrClusters=="C11"),which(mob$SnapClusters=="Gran"))]="C11Gran"
IntClusters[intersect(which(mob$ArchrClusters=="C14"),which(mob$SnapClusters=="IN4"))]="C14IN4"
IntClusters[IntClusters %in% c("C5","C8","C2","C1","C6","C7","C3","C4")] = "else"
IntClusters[IntClusters %in% c("C10","C13","C19","C20")] = "else"
IntClusters[IntClusters %in% c("C15", "C18","C17","C12","C21","C9","C11","C14")] = "else"
IntClusters.used = unique(IntClusters)[unique(IntClusters) != "else"]
mob$IntClusters = IntClusters
p = plotEmbedding(mob, colorBy = "cellColData", name = "IntClusters", embedding = "UMAP")
plotPDF(p, name = "Plot-UMAP-IntClusters.pdf", ArchRProj = mob, addDOC = FALSE, width = 6, height = 6)

markersGS2 = getMarkerFeatures(mob, groupBy = "IntClusters", useGroups = IntClusters.used, useMatrix = "GeneScoreMatrix")
saveRDS(markersGS2, file = "markersGS.Int.rds")

#gs = t(assay(markersGS.Int)[which(elementMetadata(markersGS.Int)$name == "Dlx6"),])[,1]
#names(gs) = match(names(gs), traj.maps)
#gs = gs[order(gs, decreasing = F)]
#xs = factor(names(gs), levels = names(gs), ordered = T)
#qplot(xs, gs, size = 2, color = "red")

## Plot Slc32a1

plotPDF(
    plotMarkerHeatmap(
        seMarker = markersGS2[match(markers.N, rowData(markersGS2)$name),], 
        #cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
        labelMarkers = "Slc32a1",
        transpose = TRUE
    ),
    name = "MarkerHeatmap2-Slc32a1",
    width = 8, height = 6, ArchRProj = mob, addDOC = FALSE
)


traj.maps = strsplit('C9IN3,C11Gran,C12IN5,C12IN6,C14IN4,C15IN1,C15IN2,C15IN7,C16,C17Gran,C17IN1,C17IN2,C18Gran,C18IN2,C21IN8',',')[[1]]

tryTraj = function(trajs){
    trajectories = lapply(trajs, function(s) strsplit(s, '_')[[1]])
    names(trajectories) = trajs
    trajectories = lapply(trajectories, function(t){traj.maps[as.integer(t)]})
    for (s in names(trajectories)) {
        trajectory = trajectories[[s]]
        name = paste0('t', s)
        mob = addTrajectory(
            ArchRProj = mob, 
            name = name,
            groupBy = "IntClusters",
            trajectory = trajectory, 
            embedding = "UMAP", 
            force = TRUE
        )
        plotPDF(
            plotTrajectory(mob, trajectory = name, colorBy = "cellColData", name = name), 
            name = paste0("Trajectory-v3-", name, ".pdf"), 
            ArchRProj = mob, addDOC = FALSE, width = 5, height = 5
        )
    }
    return(mob)
}

mob = tryTraj("12_13_2")
#mob = addTrajectory(mob, "t12_14_13_2", traj.maps[c(12,14,13)], "IntClusters", embedding = "UMAP", force = TRUE)
#gsm.sub = gsm[which(elementMetadata(gsm)$name %in% c("Dcx","Crhr1","Dlx6",'Calb1','Dlx1','Gabra2','Syt6','Pcp4')),]
#saveRDS(gsm.sub, file = 'gsm.sub.rds')
#traj.cells = mob$t12_14_13
#gsm = getMatrixFromProject(mob, useMatrix = "GeneScoreMatrix")
#ts = mob$t12_14_13[!is.na(mob$t12_14_13)]
#gs = assay(gsm)[which(elementMetadata(gsm)$name == "Dcx"),][!is.na(mob$t12_14_13)]
#ts = ts[gs > 0]
#gs = gs[gs > 0]
#time_gene.df = data.frame(t=ts,g=gs)
#saveRDS(time_gene.df, "time_gene.df.rds")
#qplot(time_gene.df$t, time_gene.df$g,  color = time_gene.df$g)
trajGSM  = getTrajectory(mob, name = paste0('t','12_13_2'), useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
saveRDS(trajGSM, paste0('trajGSM-t','12_13_2','.rds'))

trajGSM <- readRDS("/path/to/trajGSM-t12_13_2.rds")
#ggdf = data.frame(x=1:100,y=assay(trajGSM)[elementMetadata(trajGSM)$name == "Crhr1",])
#ggplot(ggdf, aes(x,y))+geom_point()+geom_smooth(aes(x,y),se=F,method = 'loess')

p = plotTrajectoryHeatmap(trajGSM , varCutOff = 0.9, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p, name = "PseudoTimeHeatmaps-v3", ArchRProj = mob, addDOC = FALSE, width = 6, height = 8)
#genes = sapply(rownames(p@ht_list$GeneScoreMatrix@matrix), function(i)strsplit(i,':')[[1]][2])

#dist.r = dist(assay(trajGSM), method="euclidean")
#hc = hclust(dist.r)
#hc.clusters = cutree(hc, k = 4)
#names(hc.clusters) = sapply(names(hc.clusters), function(i){strsplit(i,':')[[1]][2]})
#saveRDS(hc.clusters, file = "hc.clusters.rds")

## Split directly
genes = sapply(rownames(p@ht_list$GeneScoreMatrix@matrix), function(i)strsplit(i,':')[[1]][2])
landmarkers = c('Prkg1', 'Gfod1', 'Sorcs3')
landmarkers.idx.start = c(1, match(landmarkers, genes))
landmarkers.idx.end = c(match(landmarkers, genes)-1,length(genes))
genes.ls = lapply(1:4, function(i){genes[seq(landmarkers.idx.start[[i]], landmarkers.idx.end[[i]])]})
names(genes.ls) = 1:4

## GO
library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)

entrezid.ls = lapply(genes.ls, function(gs)as.vector(bitr(gs, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID, "character"))
enrichgo.ls = run_enrichGO(entrezid.ls, num.cores = 1)
saveRDS(enrichgo.ls, file = "enrichGO.ls.rds")

godf = Reduce(rbind, lapply(enrichgo.ls, as.data.frame))
group = data.frame(Group = rep(1:4, sapply(enrichgo.ls, nrow)))
godf = cbind(group, godf)
write.table(godf, "godf.csv", row.names = F, sep = ',')

enrichRes.dir = "GO_tables"
if(!dir.exists(enrichRes.dir)) dir.create(enrichRes.dir)
.=lapply(names(enrichgo.ls), function(s){
    fp = sprintf("%s/GO-enrich_cluster%s.csv",enrichRes.dir,s)
    write.table(enrichgo.ls[[s]], fp, sep = ',', row.names = FALSE, quote = F)
    message(paste("[INFO] Save to", fp))
})

godf = go.df[!is.na(go.df$ID),]
write.table(godf, "go.df.csv", row.names = F, sep = ',')
