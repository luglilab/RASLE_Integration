library(DropletUtils)
library(scater)
library(scran)
library(mbkmeans)
library(HDF5Array)
# import sparse matrix
sceSLE <- read10xCounts(samples = "/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/SLE_10X/")
# import metadata
metadataSLE <- read.delim("/Users/simonepuccio/Documents/HumanitasProjects/SP023/SLE_dataset/ResultFiles/RNA_sequencing_result/SDY997_EXP15176_celseq_meta.tsv.725704")
# compute size factor
sizeFactors(sceSLE) <- colSums(assay(sceSLE))
# add metadata 
sceSLE$disease <- metadataSLE$disease
#cpm(sceSLE) <- calculateCPM(sceSLE)
#tpm(sceSLE) <- calculateCPM(sceSLE)
# get Mito and Ribo geneID 
mito_genes <- grep("^MT-",rowData(sceSLE)$Symbol)
ribo_genes <- grep("^RP[SL]",rowData(sceSLE)$Symbol)
head(mito_genes,10)
# compute QC metrics
sceSLE <- calculateQCMetrics(sceSLE, feature_controls = list(Mito = mito_genes, Ribo = ribo_genes))
# determine outlier 
qc.lib <- isOutlier(sceSLE$total_counts, log=TRUE, nmads=3, type="lower",batch=sceSLE$disease)
qc.nexprs <- isOutlier(sceSLE$total_features_by_counts, nmads=3, log=TRUE, type="lower",batch=sceSLE$disease)
qc.mito <- isOutlier(sceSLE$pct_counts_Mito, nmads=3, type="higher",batch=sceSLE$disease)
qc.ribo <- isOutlier(sceSLE$pct_counts_Ribo, nmads=3, type="higher",batch=sceSLE$disease)
#
discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
table(discard)
sceSLE$discard <- discard
#
counts <- assay(sceSLE, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sceSLE) <- log2(t(t(counts)/size.factors) + 1)
assayNames(sceSLE)
# plot
plotColData(sceSLE, x = "disease", y="total_counts", colour_by = "discard") + scale_y_log10() + ggtitle("Total count")
plotColData(sceSLE, x = "disease", y="total_features_by_counts", colour_by = "discard") + scale_y_log10() + ggtitle("Detected features")
plotColData(sceSLE, x = "disease", y="pct_counts_Mito", colour_by = "discard") + ggtitle("Mito Percent")
plotColData(sceSLE, x = "disease", y="pct_counts_Ribo", colour_by = "discard") + ggtitle("Ribo Percent")
# discard 
# remove 
sceSLE <- sceSLE[,!discard]
sceSLE
# get number of gene that row sum is 0
table(rowSums(counts(sceSLE)) == 0)
# num_umis <- 1
# num_cells <- 0.05*ncol(sceSLE)
# is_expressed <- rowSums(counts(sceSLE) >= num_umis ) >= num_cells
# sceSLE <- sceSLE[is_expressed,]
# sceSLE
# compute libr size factor 
lib.sf <- librarySizeFactors(sceSLE)
summary(lib.sf)
hist(log10(lib.sf), xlab="Log10 Size factor", col='grey80')
#  mini-batch k-means algorithm
set.seed(243)
mbkm <- mbkmeans(sceSLE, clusters = 10, reduceMethod = NA)
table(mbkm$Clusters)
# compute sum factor
scran.sf <- computeSumFactors(sceSLE, cluster=mbkm$Clusters,sf.out=T)
plot(lib.sf, scran.sf, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=mbkm$Clusters)
abline(a=0, b=1, col="red")
# get row counts
counts <- assay(sceSLE, "counts")
# get log counts
logcounts(sceSLE) <- log2(t(t(counts)/lib.sf) + 1)
# print 
sceSLE
logcounts(sceSLE)
#
set.seed(333)
dbl.out <- doubletCluster(sceSLE, mbkm$Clusters)
rownames(dbl.out)[isOutlier(dbl.out$N, type="lower", nmads=3, log=TRUE)]
# save 
saveHDF5SummarizedExperiment(sceSLE, dir="/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/RA_10X/sce_SLE",
                             prefix = "sceSLE",
                             replace = TRUE)
# export workspace 
save.image(file='/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/SLEqualitycontrol.RData')
