library(SnapATAC)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"
xsp = readRDS(paste0(PROJECT_DIR, "/xsp.rds"))

## Or read from file `cluster_anno.csv`
cinfo = readRDS(paste0(PROJECT_DIR, "/cinfo.rds"))

cell.type = list(
    OEC = "olfactory ensheathing cell",
    End = "endothelial cell",
    EX = "excitatory neuron",
    AST = "astrocyte",
    OLI = "oligodendrocyte",
    MG = "microglial cell",
    OSN = "olfactory sensory neuron",
    OPC = "oligodendrocyte precursor cell",
    Purk = "Purkinje cell",
    IN = "inhibitory neuron",
    Gran = "granular cell"
)

## Cell annotated table
df = readRDS(paste0(PROJECT_DIR, "/xsp__metaData.rds"))
df$cluster = as.character(readRDS(paste0(PROJECT_DIR, "/xsp__cluster.rds")))
df$cell.type = cinfo$cell.type[match(df$cluster, cinfo$ordered.num)]
df$full.cell.type = as.character(cell.type[df$cell.type])
write.csv(df, "cell_anno.csv", quote = F, row.names = F)
saveRDS(df, "cell_anno.rds")