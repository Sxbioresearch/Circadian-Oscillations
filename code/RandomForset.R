###Use phenotypes predict rhythm patterns based on RandomForests
library(UBL)
library(dplyr)
library(randomForest)
library(caret)
library(rfPermute)
library(ggplot2)
library(splines)
library(SiZer)

##Use smote method to fix the data imbalance
OP_phenotype_RandomForest <- read.table("phenotype_raw.csv",header=T,sep = ",",stringsAsFactors = FALSE)
OP_phenotype_RandomForest$group <- as.factor(OP_phenotype_RandomForest$group)
for (i in 2:ncol(OP_phenotype_RandomForest)) {
  OP_phenotype_RandomForest[,i] <- as.factor(OP_phenotype_RandomForest[,i])
}
table(OP_phenotype_RandomForest$group)
set.seed(123)
OP_phenotype_RandomForest_UBL <- SmoteClassif(group ~ .,OP_phenotype_RandomForest,C.perc = list(anrhythmicity = 4.5,cluster1 = 24,cluster2=80,mix=50),k = 5,dist = "Overlap")
table(OP_phenotype_RandomForest_UBL$group)
##Use simulate data as train data,raw data as test data
train <- anti_join(OP_phenotype_RandomForest_UBL, OP_phenotype_RandomForest)
test <- semi_join(OP_phenotype_RandomForest_UBL, OP_phenotype_RandomForest)
table(train$group)
table(test$group)
##RandomForest
##mtry selection
features <- train[, -1]
target <- train$group
num_trees <- 500
ctrl <- trainControl(method = "cv", number = 5) 
param_grid <- expand.grid(mtry = seq(1, 67, by = 1))
set.seed(123)
rf_cv <- train(features, target, method = "rf", trControl = ctrl, tuneGrid = param_grid)
print(rf_cv)
##best mtry = 16
##RandomForest with significance levels of feature importance
set.seed(123)
otu_rfP_all <- rfPermute(group~., data = train, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1,mtry=16)
##result of feature importance
importance_all <- as.data.frame(importance(otu_rfP_all,scale = T))
write.csv(importance_all,"smote_importance_all.csv",quote = F)
##key phenotype select
set.seed(123)
##cross validation
OP_all_both.cv <- replicate(5,rfcv(train[-1], train$group,cv.fold = 10,mtry=identity,scale="new",step=-1), simplify = FALSE)
OP_all_both.cv <- data.frame(sapply(OP_all_both.cv, '[[', 'error.cv'))
OP_all_both.cv$phenotype <- rownames(OP_all_both.cv)
OP_all_both.cv <- reshape2::melt(OP_all_both.cv, id = 'phenotype')
OP_all_both.cv$phenotype <- as.numeric(as.character(OP_all_both.cv$phenotype))
write.csv(OP_all_both.cv,"OP_smote_step_select.csv",quote = F)
##turning point
model <- piecewise.linear(x=OP_all_both.cv$phenotype,y=OP_all_both.cv$value,CI=TRUE,bootstrap.samples=1000,sig.level=0.05)
model
plot(model,xlab='Phenotype',ylab='Cross-validation error')


