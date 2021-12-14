Sys.setenv(LANG="en")
library(RaceID)
#require(destiny)
#require(FateID)
require(Matrix)


print("scRNA1 smartseq")
#B1 <- read.csv("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",sep="\t", header=TRUE, row.names=1)
#print(B1[1:5,1:5])
#B1 <- as.matrix(B1)
cat("nrow = ", nrow(B1), "ncol = ", ncol(B1), "\n");

print("scRNA2 10x")
#B2 <- read.csv("new10X.tsv",sep="\t", header=TRUE,row.names=1)
#B2 <- as.matrix(B2)
cat("nrow = ", nrow(B2), "ncol = ", ncol(B2), "\n");

print("MERFISH 1")
#M1 <- read.csv("MERFISH_gene_by_cell_counts1.csv",sep=",", header=TRUE,row.names=1)
cat("nrow = ", nrow(M1), "ncol = ", ncol(M1), "\n");

print("MERFISH 1")
M2 <- read.csv("newMERFISH2.csv",sep=",", header=TRUE,row.names=1)
cat("nrow = ", nrow(M2), "ncol = ", ncol(M2), "\n");




nB1 = as.data.frame(as.matrix(B1))
nB2 <- as.matrix(B2)
nB3 <- as.matrix(M1)
nB4 <- as.matrix(M2)
class(nB1)
class(nB2)
class(nB3)
class(nB4)


gene1 <- rownames(B1) 
gene2 <- rownames(B2) 
gene3 <- rownames(M1) 
gene4 <- rownames(M2) 


cat("nrow = B1 ", nrow(B1), "ncol = ", ncol(B1), "\n");
cat("nrow = B2 ", nrow(B2), "ncol = ", ncol(B2), "\n");#%identical(rownames(B1), rownames(B2))
com1 <- intersect(gene1, gene2)
com2 <- intersect(gene3, gene4)
common_genes <- intersect(com1, com2)
#Rewrite B1 and B2 with the common genes


pB1 <- B1[common_genes,]
pB2 <- B2[common_genes,]
pB3 <- M1[common_genes,]
pB4 <- M2[common_genes,]

#common_genes

x <- cbind(pB1,pB2)
cat("x nrow= ", nrow(x), "ncol = ", ncol(x), "\n");

#identical(rownames(nB1), rownames(nB2))
dim(pB1)
dim(pB2)


D <- list()
D[["singlecell1"]] <- pB1[grep("xyz",rownames(pB1),invert=TRUE),-1]
D[["singlecell2"]] <- pB2[grep("xyz",rownames(pB2),invert=TRUE),-1]
D[["seqFISH1"]] <- pB3[grep("xyz",rownames(pB3),invert=TRUE),-1]
D[["seqFISH2"]] <- pB4[grep("xyz",rownames(pB4),invert=TRUE),-1]
d0 <- cbind(D["singlecell1"],D["singlecell2"])
dY <- cbind(D["seqFISH1"],D["seqFISH2"])
d0Y <- cbind(D[["singlecell1"]],D[["singlecell2"]],D[["seqFISH1"]],D[["seqFISH2"]])

dim(d0Y)
d0Y[1:2,1:2]



sc <- SCseq(d0Y)
sc <- filterdata(sc,minexpr = 5,mintotal=5, CGenes=grep("^(mt|Gm\\d|Rp(l|s))",rownames(d0Y),value=TRUE)) #

expData <- getExpData(sc)
#types <- sub("[B0ACTGN]+\\_","",colnames(expData))
#names(types) <- colnames(expData)
dim(expData)

types <- colnames(expData)
types[grepl("MGH",types)]<- "batch1"
types[grepl("BT",types)]<- "batch1"
types[grepl("p",types)]<- "batch2"
types[grepl("cell",types)]<- "batch3"
types[grepl("q",types)]<- "batch4"

names(types) <- colnames(expData)

length(unique(types))



res <- pruneKnn(expData,knn=60,no_cores=4,batch=types,bmethod="harmony")

cl <- graphCluster(res,pvalue=0.01,min.size=5,leiden.resolution=1.0,use.leiden=TRUE)

probes <- transitionProbs(res,cl)
sc <- updateSC(sc,res=res,cl=cl)
sc <- comptsne(sc,perplexity=200)
sc <- compumap(sc,n_neighbors=50,min_dist=0.5)


#tsne
plotmap(sc)
#tsne with names 
plotsymbolsmap(sc,types)
#plot umap 
plotmap(sc,um=TRUE)
plotsymbolsmap(sc,types,um=TRUE)
#plot expression of a given gene 
#plotexpmap(sc,"Slamf1")


unique(cl$partition)
length(unique(cl$partition))
#head(cl$fr)
write.csv(as.data.frame(cl$partition), "batch_harmony_cluster.csv",row.names = TRUE, col.names = TRUE)







