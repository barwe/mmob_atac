library(SnapATAC)
library(leiden)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
xsp = readRDS(paste0(PROJECT_DIR, "/xsp.rds"))

## Select principal component
xsp = runDiffusionMaps(xsp, "bmat", num.eigs=50)

plotDimReductPW(obj=xsp, eigs.dims=1:50, down.sample=5000,
                point.size=0.3, point.color="grey", 
                point.shape=19, point.alpha=0.6,
                pdf.height=7, pdf.width=7, 
                pdf.file.name="dim_reduct_PW.pdf")

## Run cluster
xsp = runKNN(xsp, eigs.dims = 1:20, k = 15)

xsp = runCluster(
    obj = xsp,
    tmp.folder = tempdir(),
    louvain.lib = 'leiden',
    resolution = 0.5,
    seed.use=10
)

print(length(levels(xsp@cluster)))

## Visualize

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
    output.filename = sprintf("viz_cluster_%s.pdf", type)
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

saveRDS(xsp, file=paste0(PROJECT_DIR, "/xsp.rds"))

## Run hierarchical clustering

local({
    fp = 'hclust.pdf'
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