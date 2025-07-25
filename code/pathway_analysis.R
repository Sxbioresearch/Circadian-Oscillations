### Enrichment analysis
library(tidyr)

##Pathway hierarchy information
Metacyc <- read.table("MetaCyc_pathway_map.tsv",header = T,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
names(Metacyc)[1] <- "name1"
Metacyc$name1 <- gsub("\"", "", Metacyc$name1)
OP_pathway <- read.table("pathway_all.tsv",header = T,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
OP_pathway$pathway_raw <- OP_pathway$pathway
OP_pathway$pathway <- gsub("\"", "", OP_pathway$pathway)
OP_pathway <- separate(OP_pathway, col = 1, into = c("name1", "name2"), sep = ": ", extra = "merge")
OP_name1 <- OP_pathway$name1
OP_class <- subset(Metacyc, Metacyc$name1 %in% OP_name1)
OP_merge_class <- merge(OP_pathway,OP_class, by = "name1",all = T)
OP_rhythem <- read.table("pathway_rhythem.tsv",header = T,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
OP_rhythem <- OP_rhythem$name
OP_rhythem_class <-subset(OP_pathway, OP_pathway$pathway %in% OP_rhythem)
OP_anrhythem_class <- subset(OP_pathway, !(OP_pathway$pathway %in% OP_rhythem))

#enrichment
OP_class2 <- read.table("pathway_class.tsv",header = T,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
OP_class2$p <- 1
for (i in 1:nrow(OP_class2)) {
  OP_class2[i,9] <- phyper(OP_class2[i,8]-1, OP_class2[i,6], OP_class2[i,5]-OP_class2[i,6], OP_class2[i,7],lower.tail=FALSE)
}
OP_class2$BH <- p.adjust(OP_class2$p,method = "BH")
write.table(OP_class2,"enrichment_result.tsv",sep = "\t",quote = F,row.names = F)


###Sankey diagram(sum refer to reletive abundance)
OP_paint1 <- read.table("OP_paint1.tsv",header = T,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
library(tidyr)
OP_paint2 <-  separate(OP_paint1, col = 1, into = c("pathway", "name"), sep = ":", extra = "merge")
OP_paint2$count <- 1
OP_paint3 <- as.data.frame(cbind(OP_paint2$pathway,OP_paint2$species,OP_paint2$Superclass2,OP_paint2$ownership_new,OP_paint2$count,OP_paint2$sum))
names(OP_paint3) <- c("pathway","taxa","superclass2","label","count","sum")
###networkD3
library(networkD3)
library(dplyr)
links_new <- as.data.frame(OP_paint3[,2])
links_new <- cbind(links_new,OP_paint3[,1])
links_new <- cbind(links_new,OP_paint3$sum)
names(links_new) <- c("source","target","abundance")

nodes_new <- data.frame(
  name=c(as.character(links_new$source),
         as.character(links_new$target)) %>% unique()
)

OP_cluster <-  read.table("OP_cluster.csv",header = T,sep = ",",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
names(OP_cluster)[1] <- "source"
links_merge <- merge(links_new,OP_cluster,by="source",all.x = T)

group <- read.table("pathway_taxa_cluster.csv",header = T,sep = ",",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
nodes_merge <- merge(nodes_new,group,by="name",all.x = T)
#write.csv(nodes_merge,"nodes_merge.csv",quote = F)

nodes_merge <- read.table("nodes_merge.csv",header = T,sep = ",",stringsAsFactors = FALSE,quote = "",check.names=FALSE)
links_merge$IDsource <- match(links_merge$source, nodes_merge$name)-1
links_merge$IDtarget <- match(links_merge$target, nodes_merge$name)-1

colorfunc <- JS('colorfunc = function(i) { 
                  if (i == "cluster1") {
                    return "#725663";
                  } else if (i == "cluster2") {
                    return "#D49464";
                  } else if (i == "mix") {
                    return "#5B8FA8";
                  } else {
                    return d3.schemeCategory20[Math.floor(Math.random()*20)];
                  }
                };')
p <-  sankeyNetwork(Links = links_merge, Nodes = nodes_merge,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "abundance", 
                    NodeID = "name",
                    NodeGroup="group",LinkGroup="ownership_new",
                    colourScale = colorfunc,
                    iterations=0, 
                    sinksRight=FALSE)
saveNetwork(p,"OP_taxa_pathway_networkD3_all.html")
#install.packages("webshot")
library(webshot)
if(!is_phantomjs_installed()){
  install_phantomjs()
}
is_phantomjs_installed()
## [1] TRUE
webshot("OP_taxa_pathway_networkD3_all.html" , "OP_taxa_pathway_networkD3_all.pdf",vwidth = 992,
        vheight = 1000)