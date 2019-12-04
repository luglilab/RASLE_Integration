library(DropletUtils)
library(scater)
library(scran)
library(mbkmeans)
library(HDF5Array)
# import sparse matrix
sceRA <- read10xCounts(samples = "/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/RA_10X")
# import metadata
metadataRA <- read.delim("~/Documents/HumanitasProjects/SP023/RA_dataset/ResultFiles/RNA_sequencing_result/celseq_meta.tsv.725591")
# compute size factor
sizeFactors(sceRA) <- colSums(assay(sceRA))
# add metadata 
sceRA$disease <- metadataRA$disease
#cpm(sceRA) <- calculateCPM(sceRA)
#tpm(sceRA) <- calculateCPM(sceRA)
# get Mito and Ribo geneID 
mito_genes <- grep("^MT-",rowData(sceRA)$Symbol)
ribo_genes <- grep("^RP[SL]",rowData(sceRA)$Symbol)
head(mito_genes,10)
# compute QC metrics
sceRA <- calculateQCMetrics(sceRA, feature_controls = list(Mito = mito_genes, Ribo = ribo_genes))
# determine outlier 
qc.lib <- isOutlier(sceRA$total_counts, log=TRUE, nmads=3, type="lower",batch=sceRA$disease)
qc.nexprs <- isOutlier(sceRA$total_features_by_counts, nmads=3, log=TRUE, type="lower",batch=sceRA$disease)
qc.mito <- isOutlier(sceRA$pct_counts_Mito, nmads=3, type="higher",batch=sceRA$disease)
qc.ribo <- isOutlier(sceRA$pct_counts_Ribo, nmads=3, type="higher",batch=sceRA$disease)
# plot
plotColData(sceRA, x = "disease", y="total_counts", colour_by = "discard") + scale_y_log10() + ggtitle("Total count")
plotColData(sceRA, x = "disease", y="total_features_by_counts", colour_by = "discard") + scale_y_log10() + ggtitle("Detected features")
plotColData(sceRA, x = "disease", y="pct_counts_Mito", colour_by = "discard") + ggtitle("Mito Percent")
plotColData(sceRA, x = "disease", y="pct_counts_Ribo", colour_by = "discard") + ggtitle("Ribo Percent")
# discard 
discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
table(discard)
# save discarded 
sceRA$discard <- discard
# remove 
sceRA <- sceRA[,!discard]
sceRA
# get number of gene that row sum is 0
table(rowSums(counts(sceRA)) == 0)
# num_umis <- 1
# num_cells <- 0.05*ncol(sceRA)
# is_expressed <- rowSums(counts(sceRA) >= num_umis ) >= num_cells
# sceRA <- sceRA[is_expressed,]
# sceRA
# compute libr size factor 
lib.sf <- librarySizeFactors(sceRA)
summary(lib.sf)
hist(log10(lib.sf), xlab="Log10 Size factor", col='grey80')
#  mini-batch k-means algorithm
set.seed(243)
mbkm <- mbkmeans(sceRA, clusters = 10, reduceMethod = NA)
table(mbkm$Clusters)
# compute sum factor
scran.sf <- computeSumFactors(sceRA, cluster=mbkm$Clusters,sf.out=T)
plot(lib.sf, scran.sf, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=mbkm$Clusters)
abline(a=0, b=1, col="red")
# get row counts
counts <- assay(sceRA, "counts")
# get log counts
logcounts(sceRA) <- log2(t(t(counts)/lib.sf) + 1)
# print 
sceRA
logcounts(sceRA)
#
set.seed(333)
dbl.out <- doubletCluster(sceRA, mbkm$Clusters)
rownames(dbl.out)[isOutlier(dbl.out$N, type="lower", nmads=3, log=TRUE)]
# save 
saveHDF5SummarizedExperiment(sceRA, dir="/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/RA_10X/sce_RA",
                             prefix = "sceRA",
                             replace = TRUE)
# export workspace 
save.image(file='/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/RAqualitycontrol.RData')
