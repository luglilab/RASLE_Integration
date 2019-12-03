library(Matrix)
library(DropletUtils)
# import matrix
umi_counts <- read.delim("/Users/simonepuccio/Documents/HumanitasProjects/SP023/SLE_dataset//ResultFiles/RNA_sequencing_result/umi_counts.txt",
                         row.names = 1)
# prepare sparse matrix
sparse_umi_counts <- Matrix(as.matrix(umi_counts),sparse = T)
# extrapulate barcode 
barcodes <- colnames(sparse_umi_counts)
# generate Counts in the format produced by the CellRanger 
write10xCounts("/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/SLE_10X",
               sparse_umi_counts,
               barcodes =barcodes,
               gene.id = row.names(sparse_umi_counts),
               version = "3")
# save space 
rm(umi_counts)
# workspace
save.image(file='/Users/simonepuccio/Documents/HumanitasProjects/SP023/Reanalysis_RA_SLE/SLEconversion.RData')