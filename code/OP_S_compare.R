###Comparison of the rhythem species between two sites
library(UpSetR)
library(VennDiagram) 
all_merge_upset <- read.table("all_merge_upset_new_version2.csv",header = T,sep = ",",stringsAsFactors = FALSE,quote = "",check.names=FALSE,na.strings="")
all_merge_upset_list <- list(na.omit(all_merge_upset[,1]), na.omit(all_merge_upset[,2]), na.omit(all_merge_upset[,3]), na.omit(all_merge_upset[,4]), na.omit(all_merge_upset[,5]), na.omit(all_merge_upset[,6]),na.omit(all_merge_upset[,7]),
                             na.omit(all_merge_upset[,8]))  
names(all_merge_upset_list) <- colnames(all_merge_upset[1:8])    
upset(fromList(all_merge_upset_list),  
      nsets = 100,     
      nintersects = 40, 
      order.by = "freq", 
      keep.order = F, 
      mb.ratio = c(0.6,0.4),  
      #  text.scale = 2 ,
      #     queries = list(
      #list(
      # query = intersects, 
      #params = list("OP_rhythem", "S_rhythem","OP_exist","S_exist"), 
      #color = "#725663", 
      #active = T)
      #  )
)


###Comparison of the rhythem species between two sites in every individual
###85 high-frequency species
library(UpSetR)
library(VennDiagram) 
rhythempattern_merge_upset <- read.table("rhythempattern_merge_upset_new_version2.csv",header = T,sep = ",",stringsAsFactors = FALSE,quote = "",check.names=FALSE,na.strings="")

common_elements <- intersect(rhythempattern_merge_upset$OP_absence, rhythempattern_merge_upset$S_absence)

rhythempattern_merge_upset_new <- rhythempattern_merge_upset

rhythempattern_merge_upset_new1 <- rhythempattern_merge_upset_new$OP_absence[!rhythempattern_merge_upset_new$OP_absence %in% common_elements]

rhythempattern_merge_upset_new2 <- rhythempattern_merge_upset_new$S_absence[!rhythempattern_merge_upset_new$S_absence %in% common_elements]


rhythempattern_merge_upset_list_new <- list(na.omit(rhythempattern_merge_upset_new1), na.omit(rhythempattern_merge_upset_new[,2]), na.omit(rhythempattern_merge_upset_new[,3]), na.omit(rhythempattern_merge_upset_new2), na.omit(rhythempattern_merge_upset_new[,5]), na.omit(rhythempattern_merge_upset_new[,6]),
                                            na.omit(rhythempattern_merge_upset_new[,7]),
                                            na.omit(rhythempattern_merge_upset_new[,8])
                                            #na.omit(rhythempattern_merge_upset[,9]),
                                            #na.omit(rhythempattern_merge_upset[,10])
)  
names(rhythempattern_merge_upset_list_new) <- colnames(rhythempattern_merge_upset_new[1:8])    

upset(fromList(rhythempattern_merge_upset_list_new),  
      nsets = 100,     
      nintersects = 40, 
      order.by = "freq", 
      keep.order = F, 
      mb.ratio = c(0.6,0.4),  
      queries = list(
        list(
          query = intersects, 
          params = list("pattern_diff","OP_rhythem","S_rhythem","OP_exist","S_exist"), 
          color = "#D49464", 
          active = T),
        list(
          query = intersects, 
          params = list("pattern_same","OP_rhythem","S_rhythem","OP_exist","S_exist"), 
          color = "#725663", 
          active = T),
        list(
          query = intersects, 
          params = list("OP_rhythem","OP_exist","S_exist"), 
          color = "#ADB17D", 
          active = T),
        list(
          query = intersects, 
          params = list("S_rhythem","OP_exist","S_exist"), 
          color = "#5B8FA8", 
          active = T)
      )
)