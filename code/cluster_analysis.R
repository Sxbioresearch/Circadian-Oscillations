###Clustering of rhythm species

library(Mfuzz)
library(tibble)
library(geocmeans)

# The first round of clustering of oscillatory patterns
OP_target_use <- read.table("OP_rhythm_data.csv",header = T,sep = ",",stringsAsFactors = FALSE)
OP_merge_use1 <- cbind(OP_target_use[,1],OP_target_use[,4],OP_target_use[,5:12])
names(OP_merge_use1)[1:2] <- c("taxa","people")
OP <- tidyr::unite(OP_merge_use1, "name", taxa, people,remove = T)
row.names(OP) <- NULL
OP <- OP %>% column_to_rownames("name")
OP <- as.matrix(OP)
mfuzz_OP <- new('ExpressionSet',exprs = OP)
#Preprocess missing or outlier values
mfuzz_OP <- filter.NA(mfuzz_OP, thres = 0.25)
mfuzz_OP <- fill.NA(mfuzz_OP, mode = 'mean')
mfuzz_OP <- filter.std(mfuzz_OP, min.std = 0)
mfuzz_OP <- standardise(mfuzz_OP)
#cluster number selection
cluster_num <- NULL
mfuzz_cluster_OP <- NULL
for (i in 2:10) {
  set.seed(123)
  cluster_num[i] <- i
  mfuzz_cluster_OP[[i]] <- mfuzz(mfuzz_OP, c = cluster_num[i], m = mestimate(mfuzz_OP))
}
calindex <- NULL
calindex_all <- NULL
for (i in 2:10) {
  calindex[[i]] <- as.data.frame(calcqualityIndexes(OP_target_use,mfuzz_cluster_OP[[i]]$membership, m=mestimate(mfuzz_OP),indices = c("Silhouette.index", "Partition.entropy", "Partition.coeff","XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia","DaviesBoulin.index", "CalinskiHarabasz.index","Negentropy.index")))
  temp <- calindex[[i]]
  calindex_all <- rbind(calindex_all,temp)
}

#cluster analysis
set.seed(123)
cluster_num <- 2
mfuzz_cluster_OP <- mfuzz(mfuzz_OP, c = cluster_num, m = mestimate(mfuzz_OP))

#plot
#cluster1
x <- c(9,13,17,21,33,37,41,45)
y <- mfuzz_cluster_OP[[1]][1,]
mfuzz.plot2(mfuzz_OP, cl = mfuzz_cluster_OP, mfrow = c(2, 3.5), colo = "fancy",time.labels = colnames(OP),time.points=c(9,13,17,21,33,37,41,45),single = 1)
#lines(mfuzz_cluster_OP[[1]][1,],col=centre.col,lwd=centre.lwd)
lines(x, y, col = "black", lty = 1,lwd = 2)
abline(v = 24, col = "black", lty = 2)
#cluster2
y <- mfuzz_cluster_OP[[1]][2,]
mfuzz.plot2(mfuzz_OP, cl = mfuzz_cluster_OP, mfrow = c(2, 3.5), colo = "fancy",time.labels = colnames(OP),time.points=c(9,13,17,21,33,37,41,45),single = 2)
#lines(mfuzz_cluster_OP[[1]][1,],col=centre.col,lwd=centre.lwd)
lines(x, y, col = "black", lty = 1,lwd = 2)
abline(v = 24, col = "black", lty = 2)
#plot membership >0.8
#cluster1
x <- c(9,13,17,21,33,37,41,45)
y <- mfuzz_cluster_OP[[1]][1,]
mfuzz.plot2(mfuzz_OP, cl = mfuzz_cluster_OP, mfrow = c(2, 3.5), colo = "fancy",time.labels = colnames(OP),time.points=c(9,13,17,21,33,37,41,45),single = 1,min.mem = 0.8)
#lines(mfuzz_cluster_OP[[1]][1,],col=centre.col,lwd=centre.lwd)
lines(x, y, col = "black", lty = 1,lwd = 2)
abline(v = 24, col = "black", lty = 2)
#cluster2
y <- mfuzz_cluster_OP[[1]][2,]
mfuzz.plot2(mfuzz_OP, cl = mfuzz_cluster_OP, mfrow = c(2, 3.5), colo = "fancy",time.labels = colnames(OP),time.points=c(9,13,17,21,33,37,41,45),single = 2,min.mem = 0.8)
#lines(mfuzz_cluster_OP[[1]][1,],col=centre.col,lwd=centre.lwd)
lines(x, y, col = "black", lty = 1,lwd = 2)
abline(v = 24, col = "black", lty = 2)
#Standardized result
write.csv(mfuzz_OP@assayData$exprs,"OP_standard.csv",quote = F)
#value of membership
write.csv(mfuzz_cluster_OP$membership,"OP_membership.csv",quote = F)
##species with membership < 0.8 repeat the above process
