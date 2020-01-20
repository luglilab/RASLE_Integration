library(DropletUtils)
library(scater)
library(scran)
library(mbkmeans)
#######################################
Ptz7Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz7_Involved/",
                              sample.names = "Ptz7Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz7Involved) <- colData(Ptz7Involved)$Barcode
#######################################
Ptz7Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz7_Uninvolved/",
                                sample.names = "Ptz7Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz7Uninvolved) <- colData(Ptz7Uninvolved)$Barcode
#######################################
Ptz8Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz8_Involved/",
                              sample.names = "Ptz8Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz8Involved) <- colData(Ptz8Involved)$Barcode
#######################################
Ptz8Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz8_Uninvolved/",
                                sample.names = "Ptz8Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz8Uninvolved) <- colData(Ptz8Uninvolved)$Barcode
#######################################
#######################################
Ptz10Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz10_Involved/",
                              sample.names = "Ptz10_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz10Involved) <- colData(Ptz10Involved)$Barcode
#######################################
Ptz10Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz10_Uninvolved/",
                                sample.names = "Ptz10_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz10Uninvolved) <- colData(Ptz10Uninvolved)$Barcode
#######################################
#######################################
Ptz11Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz11_Involved/",
                              sample.names = "Ptz11_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz11Involved) <- colData(Ptz11Involved)$Barcode
#######################################
Ptz11Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz11_Uninvolved/",
                                sample.names = "Ptz11_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c( "2"), genome = NULL)
colnames(Ptz11Uninvolved) <- colData(Ptz11Uninvolved)$Barcode
#######################################
#######################################
Ptz12Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz12_Involved/",
                              sample.names = "Ptz12_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz12Involved) <- colData(Ptz12Involved)$Barcode
#######################################
Ptz12Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz12_Uninvolved/",
                                sample.names = "Ptz12_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz12Uninvolved) <- colData(Ptz12Uninvolved)$Barcode
#######################################
#######################################
Ptz13Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz13_Involved/",
                              sample.names = "Ptz13_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz13Involved) <- colData(Ptz13Involved)$Barcode
#######################################
Ptz13Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz13_Uninvolved/",
                                sample.names = "Ptz13_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz13Uninvolved) <- colData(Ptz13Uninvolved)$Barcode
#######################################
#######################################
Ptz14Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz14_Involved/",
                              sample.names = "Ptz14_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz14Involved) <- colData(Ptz14Involved)$Barcode
#######################################
Ptz14Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz14_Uninvolved/",
                                sample.names = "Ptz14_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("2"), genome = NULL)
colnames(Ptz14Uninvolved) <- colData(Ptz14Uninvolved)$Barcode
#######################################
#######################################
Ptz15Involved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz15_Involved/",
                              sample.names = "Ptz15_Involved",
                              col.names = FALSE, 
                              type = c("auto", "sparse", "HDF5"),
                              version = "2", genome = NULL)
colnames(Ptz15Involved) <- colData(Ptz15Involved)$Barcode
#######################################
Ptz15Uninvolved <- read10xCounts("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/Ptz15_Uninvolved/",
                                sample.names = "Ptz15_Uninvolved",
                                col.names = FALSE, type = c("auto", "sparse", "HDF5"),
                                version = c("auto", "2", "3"), genome = NULL)
colnames(Ptz15Uninvolved) <- colData(Ptz15Uninvolved)$Barcode
#######################################
common <- intersect(rownames(Ptz7Involved), rownames(Ptz7Uninvolved))
sce <- cbind(Ptz7Involved[common,], Ptz7Uninvolved[common,],
             Ptz8Involved[common,], Ptz8Uninvolved[common,],
             Ptz10Involved[common,], Ptz10Uninvolved[common,],
             Ptz11Involved[common,], Ptz11Uninvolved[common,],
             Ptz12Involved[common,], Ptz12Uninvolved[common,],
             Ptz13Involved[common,], Ptz13Uninvolved[common,],
             Ptz14Involved[common,], Ptz14Uninvolved[common,],
             Ptz15Involved[common,], Ptz15Uninvolved[common,])
write.table(colnames(sce),"/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/metadataIDB.txt",quote = F,
            row.names = F)
setwd("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/")
system("awk \'BEGIN{FS=OFS=\"\t\";print \"CellID\",\"Samples\",\"Disease\",\"Tissue\",\"Status\",\"PatientID\",\"Chemistry\"}{if($1 ~/Ptz7I/) print $1,\"128\",\"IBD\",
       \"ILEUM\",\"Involved\",\"Ptz7\",\"V2\" ;else if($1 ~/Ptz7U/) print $1,\"129\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz7\",\"V2\" ; else if($1 ~/Ptz8I/) print $1,\"138\",
       \"128\",\"IBD\",\"ILEUM\",\"Involved\",\"Ptz8\",\"V2\"; else if($1 ~/Ptz8U/) print $1,\"135\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz8\",\"V2\";else if($1 ~/Ptz10I/) print $1,
       \"158\",\"IBD\",\"ILEUM\",\"Involved\",\"Ptz10\",\"V2\";else if($1 ~/Ptz10U/) print $1,\"159\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz10\",\"V2\";else if($1 ~/Ptz11I/) print $1,
       \"181\",\"IBD\",\"ILEUM\",\"Involved\",\"Ptz11\",\"V2\" ;else if($1 ~/Ptz11U/) print $1,\"180\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz11\",\"V2\";else if($1 ~/Ptz12I/) print $1,
       \"187\",\"IBD\",\"ILEUM\",\"Involved\",\"Ptz12\",\"V2\";else if($1 ~/Ptz12U/) print $1,\"186\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz12\",\"V2\";else if($1 ~/Ptz12I/) print $1,
       \"187\",\"IBD\",\"ILEUM\",\"Involved\",\"Ptz12\",\"V2\" ;else if($1 ~/Ptz12U/) print $1,\"186\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz12\",\"V2\";else if($1 ~/Ptz13I/) print $1,\"190\",
       \"IBD\",\"ILEUM\",\"Involved\",\"Ptz13\",\"V2\" ;else if($1 ~/Ptz13U/) print $1,\"189\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz13\",\"V2\";else if($1 ~/Ptz14I/) print $1,\"193\",\"IBD\",
       \"ILEUM\",\"Involved\",\"Ptz14\",\"V2\";else if($1 ~/Ptz14U/) print $1,\"192\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz14\",\"V2\";else if($1 ~/Ptz15I/) print $1,\"196\",\"IBD\",\"ILEUM\",
       \"Involved\",\"Ptz15\",\"V2\" ;else if($1 ~/Ptz15U/) print $1,\"195\",\"IBD\",\"ILEUM\",\"Uninvolved\",\"Ptz8\",\"V2\"}' metadataIDB.txt > metadata_tableIDB.txt ")
metadata_tableIDB <- read.delim("/mnt/hpcserver1_datadisk2_spuccio/SP023_Erc_preliminary/IBD_Analysis/GSE134809_RAW/metadata_tableIDB.txt")
metadata(sce)$Samples <-  metadata_tableIDB$Samples
metadata(sce)$Disease <-  metadata_tableIDB$Disease
metadata(sce)$Tissue <-  metadata_tableIDB$Tissue
metadata(sce)$Status <-  metadata_tableIDB$Status
metadata(sce)$PatientID <-  metadata_tableIDB$PatientID
metadata(sce)$Chemistry <-  metadata_tableIDB$Chemistry
#####  Computing barcode ranks 
set.seed(1234)
br.out <- barcodeRanks(counts(sce))
# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
#####
rowData(sce)$ensembl_gene_id <- rownames(sce)
sce <- getBMFeatureAnnos(sce, 
                         filters = "ensembl_gene_id",
                         attributes = c("ensembl_gene_id", "hgnc_symbol", 
                                        "start_position", "end_position", "chromosome_name"),
                         dataset = "hsapiens_gene_ensembl")

mito_genes <- grep("^MT-",rowData(sce)$Symbol)
ribo_genes <- grep("^RP[SL]",rowData(sce)$Symbol)
head(mito_genes,10)
head(ribo_genes,10)
sce <- calculateQCMetrics(sce, feature_controls = list(Mito = mito_genes, Ribo = ribo_genes),BPPARAM=MulticoreParam(10))
rm(Ptz7Involved)
rm(Ptz7Uninvolved)
rm(Ptz8Involved)
rm(Ptz8Uninvolved)
rm(Ptz10Involved)
rm(Ptz10Uninvolved)
rm(Ptz11Involved)
rm(Ptz11Uninvolved)
rm(Ptz12Involved)
rm(Ptz12Uninvolved)
rm(Ptz13Involved)
rm(Ptz13Uninvolved)
rm(Ptz14Involved)
rm(Ptz14Uninvolved)
rm(Ptz15Involved)
rm(Ptz15Uninvolved)
#########################
rownames(sce) <- rowData(sce)$Symbol
qcstats <-perCellQCMetrics(sce, subsets=list(Mito=grep("MT-", rownames(sce))),BPPARAM=MulticoreParam(10))
colData(sce) <- cbind(colData(sce), qcstats)
########################
keep.total <- isOutlier(sce$sum, type="lower")
keep.total <- isOutlier(sce$sum, type="lower")
filtered <- sce[,keep.total]
plotHighestExprs(sce, exprs_values = "counts")
