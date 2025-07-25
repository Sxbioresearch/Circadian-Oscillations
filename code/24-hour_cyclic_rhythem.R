library(vegan)
library(phyloseq)

###alpha diversity
OP <- read.table("OP_species_absolute.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
#The function of alpha diversity
alpha_diversity <- function(x) {
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')   
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  result <- data.frame(Shannon, Simpson, goods_Coverage)
  result
}
OP1 <- t(OP)
alpha <- alpha_diversity(OP1)
write.csv(alpha, 'OP_alpha.csv', quote = FALSE)

### Autocorrelation of alpha diversity
OP_alpha = read.table("OP_alpha.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
#add group information
OP_group = read.table("OP_group.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
OP_all <- cbind(OP_alpha,OP_group)
OP_autoresult <- NULL
for (j in seq(1,176,8)) {
  a <- acf(OP_all[j:(j+7),1],lag.max=NULL,type = c("correlation"),na.action = na.pass,plot = F)
  a <- with(a,data.frame(lag,acf))
  OP_autoresult <- rbind(OP_autoresult,a)
}
#group information
OP_auto_paint <- cbind(OP_autoresult,OP_group)
write.csv(OP_auto_paint,"OP_alpha_auto.csv",quote = F)

###beta diversity
data <- read.table("OP_species_absolute_filter.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data_trans=t(data)
data_trans1 = otu_table(data_trans, taxa_are_rows = F)
data_trans2 = phyloseq(data_trans1)
BC=phyloseq::distance(data_trans2, method = "bray")
BC1 <- as.matrix(BC)
write.csv(BC1, 'OP_BC.csv', quote = FALSE)

#permanova test
dis <- read.table('OP_BC.csv',header = T, row.names = 1, sep = ',', stringsAsFactors = FALSE, check.names = FALSE) #读入距离矩阵
OP_group <- read.table('OP_group.csv', header = T, sep = ',', stringsAsFactors = FALSE)#读入分组文件
adonis_result_dis <- adonis2(dis~People+Time, OP_group, permutations = 999,by="margin") #根据 group$People 这一列样本分组信息进行 PERMANOVA 分析，随机置换检验 999 次
adonis_result_dis

