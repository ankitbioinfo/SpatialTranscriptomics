
Sys.setenv(LANG="en")

library(batchelor)



#B1 <- matrix(rnorm(10000), ncol=100) # Batch 1
#write.csv(B1,"mat1.csv", row.names = FALSE)

#B2 <- matrix(rnorm(10000), ncol=100) # Batch 2
#write.csv(B2,"mat2.csv", row.names = FALSE)

B1 <-read.csv("for_R_spatial.csv", sep=",", header=TRUE,row.names=1)
B2 <-read.csv("for_R_sc.csv", sep=",", header=TRUE,row.names=1)


B1 <- as.data.frame(t(B1[,-1]))
B2 <- as.data.frame(t(B2[,-1]))


out <- findMutualNN(B1, B2, k1=20)
head(out$first)
head(out$second)
write.csv(out,"R_match_pairing_common.csv", row.names = FALSE)
