suppressMessages(library(ArchR))
addArchRGenome("mm10")
addArchRThreads(threads = 16)

outputDir = "/home/xiangrong1/hw5_COHORT/project/mmob/ArchRAnalysis"
if(!dir.exists(outputDir)) { dir.create(outputDir) }
setwd(outputDir)

bamFiles = "/home/xiangrong1/hw5_COHORT/project/mmob/mob_atac.markdup.sorted.bam"
names(bamFiles) = "mob"

bamFlag = list(
    isPaired = TRUE,
    isProperPair = TRUE,
    isUnmappedQuery = FALSE,
    hasUnmappedMate = FALSE,
    isMinusStrand = TRUE,
    isMateMinusStrand = NA,
    isFirstMateRead = NA,
    isSecondMateRead = NA,
    isSecondaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE,
    isDuplicate = FALSE
)

ArrowFiles = createArrowFiles(
    inputFiles = bamFiles,
    sampleNames = names(bamFiles),
    filterTSS = 4,
    filterFrags = 1000,
    addTileMat = FALSE,
    addGeneScoreMat = TRUE,
    gsubExpression= ":.*",
    bamFlag = bamFlag
)