library(DropletUtils)
library(scater)
library(scran)
sceRA <- read10xCounts(samples = "/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/RA_10X")
sceSLE <- read10xCounts(samples = "/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/SLE_10X/")

cpm(sceRA) <- calculateCPM(sceRA)
tpm(sceRA) <- calculateCPM(sceRA)
sceRA <- normalize(sceRA)
clusters <- quickCluster(sceRA)
sce <- computeSumFactors(sce, clusters=clusters)


libsizes <- colSums(counts)




mt.genes <- rownames(sceRA)[grep("^MT-",rownames(sceRA))]

sceRA <- calculateQCMetrics(sceRA, feature_controls = list(mito = mt.genes))

plotScater(sceRA, block1 = "ident", nfeatures = 1000)

