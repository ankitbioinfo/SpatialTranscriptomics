Sys.setenv(LANG="en")
library(scRNAseq)
library(batchelor)
library(scran)
library(scater)


print("scRNA1 smartseq")
B1 <- read.csv("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",sep="\t", header=TRUE, row.names=1)
#print(B1[1:5,1:5])
#B1 <- as.matrix(B1)
cat("nrow = ", nrow(B1), "ncol = ", ncol(B1), "\n");

print("scRNA2 10x")
B2 <- read.csv("GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv",sep="\t", header=TRUE,row.names=1)
#B2 <- as.matrix(B2)
cat("nrow = ", nrow(B2), "ncol = ", ncol(B2), "\n");

print("MERFISH 1")
M1 <- read.csv("MERFISH_gene_by_cell_counts1.csv",sep=",", header=TRUE,row.names=1)
cat("nrow = ", nrow(M1), "ncol = ", ncol(M1), "\n");

print("MERFISH 1")
M2 <- read.csv("MERFISH_gene_by_cell_counts2.csv",sep=",", header=TRUE,row.names=1)
cat("nrow = ", nrow(M2), "ncol = ", ncol(M2), "\n");



gene1 <- rownames(B1) 
gene2 <- rownames(B2) 
gene3 <- rownames(M1) 
gene4 <- rownames(M2) 
identical(rownames(B1), rownames(B2))
com1 <- intersect(gene1, gene2)
com2 <- intersect(gene3, gene4)
common_genes <- intersect(com1, com2)
#Rewrite B1 and B2 with the common genes
pB1 <- B1[common_genes,]
pB2 <- B2[common_genes,]
pB3 <- M1[common_genes,]
pB4 <- M2[common_genes,]

#nB1 <- cbind(pB1,pB2)
#nB2 <- cbind(pB3,pB4)
#common_genes

#identical(rownames(nB1), rownames(nB2))
#dim(nB1)
#dim(nB2)



sce1 <- SingleCellExperiment(list(counts=pB1))
sce2 <- SingleCellExperiment(list(counts=pB2))
sce3 <- SingleCellExperiment(list(counts=pB3))
sce4 <- SingleCellExperiment(list(counts=pB4))
sizeFactors(sce1) <- 2^rnorm(ncol(sce1))
sizeFactors(sce2) <- 2^rnorm(ncol(sce2))
sizeFactors(sce3) <- 2^rnorm(ncol(sce3))
sizeFactors(sce4) <- 2^rnorm(ncol(sce4))

#normed <- applyMultiSCE(sce1, sce2,sce3,sce4, FUN=multiBatchNorm)
combined <- cbind(sce1, sce2,sce3,sce4)
batch <- rep(1:4, c(ncol(sce1), ncol(sce2),ncol(sce3),ncol(sce4)))
combined <- applyMultiSCE(combined, batch=batch, FUN=multiBatchNorm)
chosen.hvgs <- getTopHVGs(combined,n=200)





f.out <- fastMNN(combined,batch=batch,subset.row=chosen.hvgs)
#The corrected matrix in the reducedDims() contains the low-dimensional corrected 
#coordinates for all cells, which we will use in place of the PCs in our downstream 
#analyses.
dim(reducedDim(f.out, "corrected"))
str(reducedDim(f.out, "corrected"))
#rle(f.out$batch)
dim(f.out)

#A reconstructed matrix in the assays() contains the corrected expression values for 
#each gene in each cell, obtained by projecting the low-dimensional coordinates in corrected 
#back into gene expression space. We do not recommend using this for anything other than 
#visualization

assay(f.out, "reconstructed")
set.seed(105)
out <- runTSNE(f.out, dimred="corrected")
plotTSNE(out, colour_by="batch")

n.gr <- buildSNNGraph(f.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=f.out$batch)
tab.mnn
write.csv(cbind(batch,clusters.mnn),"MNN_batch_correction_cluster_output.csv",quote=FALSE,row.names=TRUE)


cor.exp <- assay(f.out)[1,]
hist(cor.exp, xlab="Corrected expression for gene 1", col="grey80") 






comb <- correctExperiments(combined,batch=batch, PARAM=NoCorrectParam())
set.seed(100)
comb <- runPCA(comb, subset_row=chosen.hvgs)
comb <- runTSNE(comb, dimred="PCA")
#comb <- runUMAP(comb, dimred="PCA")
plotTSNE(comb, colour_by="batch")



