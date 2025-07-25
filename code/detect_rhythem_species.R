###Detection process of high-frequency rhythmic species
library(tidyverse)
library(circacompare)
library(dplyr)

###circacompare_method

##simulation random dataset
data <- read.table("OP_species_absolute_filter.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data <- t(cycdata)
#Randomize the values of the eight time points for each species within each participant
#random1
random1 <- apply(data, 2, function(x) {
  x <- as.numeric(x)
  x_perm <- unlist(lapply(split(x, (seq_along(x)-1) %/% 8), function(y) {
    y[sample(length(y))]
  }))
  x_perm
})
random1 <- as.data.frame(random1)
rownames(random1) <- rownames(data)
write.csv(random1,"random1.csv",quote = F)
###repeat three times to obtain three random datasets

###detect rhythem
data <- read.table("OP_species_absolute_filter.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
OP_filter_trans  <- t(data)
OP_group <- read.table("OP_group_merge.csv",header = T,sep = ",",stringsAsFactors = FALSE)
OP_merge_filter <- cbind(OP_group,OP_trans_filter)
test = OP_merge_filter  %>% group_split(People)
people <- NULL
result <-  NULL
a <- NULL
for (i in 1:22) {
  people[[i]] <- test[[i]][,8:281]
  people[[i]] <- people[[i]][,colSums(people[[i]]^2) !=0]
  people[[i]] <- cbind(test[[i]][,1:7],people[[i]])
  result[[i]] <- list()
  for (j in 8:ncol(people[[i]])) {
    out <- circa_single(x = people[[i]], col_time = "Timecompare",     col_outcome = colnames(people[[i]][j]),timeout_n = 100000)
    result[[i]] <- as.data.frame(cbind(result[[i]],out[[2]]$value))
  }
  a[[i]] <- list()
  for (j in 8:ncol(people[[i]])) {
    a[[i]] <-  rbind(a[[i]],colnames(people[[i]][j]))
  }
  names(result[[i]]) <- a[[i]]
  result[[i]] <- as.data.frame(cbind(out[[2]]$parameter,result[[i]]))
  names(result[[i]])[1] <- ""
  result[[i]] <- t(result[[i]])
  for (k in 2:nrow(result[[i]])) {
    b <- as.data.frame(test[[i]][2,2])
    c <- b[1,1]
    result[[i]][k,6] <- c
  }
  result[[i]][1,6] <- "group"
}
name=paste("circa_result_",c(1:22),".csv",sep="")
for(i in 1:22){
  write.table(result[[i]],name[i],quote=FALSE,sep=',',col.names = F)
}

###empJTK_method
##simulation random dataset
#100 repeat random datasets
#set1
data <- read.table("OP_eJTK_raw_set1.txt",sep = "\t",row.names = 1,stringsAsFactors = FALSE)
names(data) <- c("17","21","33","37","41","45","57","61")
data_random <- NULL
for (i in 1:100) {
  data_random[[i]] <- as.data.frame(t(apply(data, 1, function(x) sample(x))))
  data_random[[i]]  <- tibble::rownames_to_column(data_random[[i]] , var="#")
  names(data_random[[i]]) <- c("#","17","21","33","37","41","45","57","61")
}
name=paste("OP_random_set1_",c(1:100),".txt",sep="")
for(i in 1:100){
  write.table(data_random[[i]],name[i],quote=FALSE,sep='\t',row.names = F)
}
#set2
data <- read.table("OP_eJTK_raw_set2.txt",sep = "\t",row.names = 1,stringsAsFactors = FALSE)
names(data) <- c("9","13","17","21","33","37","41","45")
library(dplyr)
data_random <- NULL
for (i in 1:100) {
  data_random[[i]] <- as.data.frame(t(apply(data, 1, function(x) sample(x))))
  data_random[[i]]  <- tibble::rownames_to_column(data_random[[i]] , var="#")
  names(data_random[[i]]) <- c("#","9","13","17","21","33","37","41","45")
}
name=paste("OP_random_set2_",c(1:100),".txt",sep="")
for(i in 1:100){
  write.table(data_random[[i]],name[i],quote=FALSE,sep='\t',row.names = F)
}
###eJTK.sh detection rhythmic species by empJTK method




