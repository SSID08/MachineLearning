#Load libraries
library(readxl)
library(mixOmics)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(caret)
library(kernlab)
library(e1071)
library(rgl)
library(tidyverse)
library(dplyr)

source('functions.R')

#Load data
sensory_score=read.csv('../Data/Sensory_score-1.csv',header = T,row.names = 1)
HPLC=read.csv('../Data/HPLC_data-1.csv',header = T,row.names = 1)
Enose=read.csv('../Data/Enose_data-1.csv',header = T,row.names = 1)
Bacterial_counts=read.csv('../Data/Bacterial_Counts-1.csv',header = T,row.names = 1)

#Create class vectors
Enose_sensory_class=sensory_score[rownames(Enose),"sensory"]
Enose_bacterial_count=Bacterial_counts[rownames(Enose),]
HPLC_sensory_class=sensory_score[rownames(HPLC),]
HPLC_bacterial_count=Bacterial_counts[rownames(HPLC),]

#Make class vectors into factors
Enose_sensory_class=as.factor(Enose_sensory_class)
HPLC_sensory_class=as.factor(HPLC_sensory_class)

#Create the PCA objects
PCA_enose=pca(Enose,ncomp = 5,scale = T)
PCA_HPLC=pca(HPLC,ncomp = 10,scale = T)

#Plot the PCA objects
plotIndiv(PCA_enose,ind.names = Enose_sensory_class,group = Enose_sensory_class,ellipse = T,col.per.group = c('Red','Blue','Green'),legend = T,title = 'PCA for Enose_data',style = 'lattice')
plotIndiv(PCA_HPLC,ind.names = HPLC_sensory_class,group = HPLC_sensory_class,ellipse = T,col.per.group = c('Red','Blue','Green'),legend = T,title='PCA for HPLC data',style = '3d')

#Merge bacterial count with sensory data
Bacterial_v_Sensory=merge(sensory_score,Bacterial_counts,by = 'row.names')%>%column_to_rownames('Row.names')
#Plot boxplot of sensory against bacterial count
ggplot(Bacterial_v_Sensory)+geom_boxplot(mapping=aes(x=sensory,y=TVC,group=sensory))+ggtitle('TVC .v Sensory')
ggplot(Bacterial_v_Sensory)+geom_boxplot(mapping=aes(x=sensory,y=Pseudomonads,group=sensory))+ggtitle('Pseudomonas v. Sensory')

#Plot of time against bacterial count
Bacterial_v_Sensory$ID=rownames(Bacterial_v_Sensory)
Bacterial_v_Sensory=Bacterial_v_Sensory%>%separate(ID,into = c('Temperature','Time'),sep='F')
Bacterial_v_Sensory=Bacterial_v_Sensory[-c(51,52),]
Bacterial_v_Sensory$Time=as.numeric(Bacterial_v_Sensory$Time)
Bacterial_v_Sensory$Temperature=factor(Bacterial_v_Sensory$Temperature,levels = c('0','5','10','15','20'))
ggplot(Bacterial_v_Sensory,mapping=aes(x=Time,y=TVC,group=Temperature,color=Temperature))+geom_line(alpha=.3)+geom_smooth(se=F,linewidth=.8,alpha=.2)
ggplot(Bacterial_v_Sensory,mapping=aes(x=Time,y=Pseudomonads,group=Temperature,color=Temperature))+geom_line(alpha=.3)+geom_smooth(se=F,linewidth=.8,alpha=.2)

#Matrix boxlpot of variables feafures in dataset

boxplot.matrix(as.matrix(Enose),main='Enose variables boxplot')
boxplot.matrix(as.matrix(HPLC),main='HPLC variables boxplot')

#Convert sensory score to factor

sensory_score$sensory=as.factor(sensory_score$sensory)

#Enose data preperation

Enose_sensory=merge_func(Enose,sensory_score)

Enose_bc_TVC=merge_func(Enose,Bacterial_counts)%>%select(-Pseudomonads)%>%rename('Bac_count'=TVC)
Enose_bc_pseudo=merge_func(Enose,Bacterial_counts)%>%select(-TVC)%>%rename('Bac_count'=Pseudomonads)


#HPLC data preperation

HPLC_sensory=merge_func(HPLC,sensory_score)

HPLC_bc_TVC=merge_func(HPLC,Bacterial_counts)%>%select(-Pseudomonads)%>%rename('Bac_count'=TVC)
HPLC_bc_pseudo=merge_func(HPLC,Bacterial_counts)%>%select(-TVC)%>%rename('Bac_count'=Pseudomonads)

#Set up tuning and training parameters
ctr=tune.control(sampling = 'cross',cross= 3,nrepeat= 10)
train_control=trainControl(method = 'repeatedcv',search = 'grid',number = 3,repeats = 10)

#Tune the KNN models

set.seed(2)
Enose_knn_classification=Knn_tune(data=Enose_sensory,k = c(1:10))
HPLC_knn_classification=Knn_tune(data=HPLC_sensory,k = c(1:10))

set.seed(234)
Enose_knn_pseudo=Knn_tune(data = Enose_bc_pseudo,k=c(1:10))
Enose_knn_Tvc=Knn_tune(data=Enose_bc_TVC,k=c(1:10))

set.seed(534)
HPLC_knn_pseudo=Knn_tune(data = HPLC_bc_pseudo,k=c(1:10))
HPLC_knn_Tvc=Knn_tune(data = HPLC_bc_TVC,k=c(1:10))


#Create Custom Random Forest Tune object

customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree",'nodesize','maxnodes'),
                                  class = rep("numeric", 4),
                                  label = c("mtry", "ntree",'nodesize','maxnodes'))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree,
               maxnodes=param$maxnodes,
               nodesize=param$nodesize)
}

#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes



#Tune RF models

Enose_rf_sensory=random_forest_tune_classification(data = Enose_sensory,mtry = c(1:8),ntree=c(100,200,500),nodesize=c(1,3,5),maxnodes = c(5,10,15),fold = 3,nrepeats = 10,seed = 1234)
HPLC_rf_sensory=random_forest_tune_classification(data=HPLC_sensory,mtry = c(2,5,10,15,19),ntree=c(100,200,500),nodesize = c(1,3,5),maxnodes = c(5,10,15),fold = 3,nrepeats = 10,seed = 123)

Enose_rf_TVC=random_forest_tune_regression(data=Enose_bc_TVC,mtry = c(1:8),ntree = c(100,200,500),nodesize = c(5,10,15),seed=450)
HPLC_rf_TVC=random_forest_tune_regression(data=HPLC_bc_TVC,mtry=c(2,5,10,15,19),ntree=c(100,200,500),nodesize = c(5,10,15),seed=500)

Enose_rf_pseudo=random_forest_tune_regression(data=Enose_bc_pseudo,mtry = c(1:8),ntree = c(100,200,500),nodesize = c(5,10,15),seed = 320)
HPLC_rf_pseudo=random_forest_tune_regression(data=HPLC_bc_pseudo,mtry = c(2,5,10,15,19),ntree=c(100,200,500),nodesize = c(5,10,15),seed = 100)

#Tune SVM models

Enose_svm_classification=svmrd_tune(data=Enose_sensory,cost = c(0.1,1,5,10,20),sigma=c(0.1,1,5,10),seed = 1234)
HPLC_svm_classification=svmrd_tune(data=HPLC_sensory,cost = c(0.1,1,5,10,20),sigma=c(0.1,1,5,10),seed=1234)

#Tune PLS-R

Enose_TVC_pslr=spls_tune(data = Enose_bc_TVC,num_lv = seq(1,8,1),ncomp = 5,seed = 23)
Enose_pseudo_pslr=spls_tune(data=Enose_bc_pseudo,num_lv = seq(1,8,1),ncomp=5,seed=345)
HPLC_TVC_pslr=spls_tune(data=HPLC_bc_TVC,num_lv = c(seq(1,18,2),18),ncomp = 8,seed = 43)
HPLC_pseudo_pslr=spls_tune(data=HPLC_bc_pseudo,num_lv = c(seq(1,18,2),18),ncomp = 8,seed=46)

#Run iterations for Enose data
Enose_100_out=iter_function(data = Enose_sensory,k = 5,mtry = 2,ntree = 500,maxnodes = 15,nodesize = 3,C = 5,sigma = .1,iter = 100,seed=424)
Enose_100_TVC=iter_function_regression(data = Enose_bc_TVC,k = 5,mtry = 2,ntree = 500,nodesize = 5,ncomp = 1,iter = 100,seed = 123)
Enose_100_pseudo=iter_function_regression(data = Enose_bc_pseudo,k = 6,mtry = 4,ntree = 200,nodesize = 5,ncomp = 1,iter = 100,seed=765)

#Run iterations for HPLC data
HPLC_100_out=iter_function(data=HPLC_sensory,k=5,mtry = 10,ntree = 200,maxnodes = 10,nodesize = 1,C = 5,sigma=.1,iter = 100,seed=257)
HPLC_100_TVC=iter_function_regression(data=HPLC_bc_TVC,k = 3,mtry = 2,ntree = 100,nodesize = 5,ncomp = 1,iter = 100,seed = 343)
HPLC_100_pseudo=iter_function_regression(data=HPLC_bc_pseudo,k=3,mtry = 5,ntree = 200,nodesize = 5,ncomp = 1,iter=100,seed=89)

#Extract summary metrics for the 100 iterations of the regression models
Enose_TVC_regression_summary=extract_metrics(metric = Enose_100_TVC)
Enose_Pseudo_regression_summary=extract_metrics(Enose_100_pseudo)
HPLC_TVC_regression_summary=extract_metrics(HPLC_100_TVC)
HPLC_pseudo_regression_summary=extract_metrics(HPLC_100_pseudo)

#Extract missclassified summary statistics for the 100 iterations of the Classification models

Enose_misclassified_df=data.frame(t(matrix(c(Enose_100_out$`Missclassified Samples KNN`,Enose_100_out$`Missclasssified samples SVM`,Enose_100_out$`Missclassified samples RF`),nrow = 3,byrow = T))/Enose_100_out$`Sum of Class samples`)
colnames(Enose_misclassified_df)=c('KNN','SVM','RF')
rownames(Enose_misclassified_df)=c('1','2','3')

HPLC_misclassified_df=data.frame(t(matrix(c(HPLC_100_out$`Missclassified Samples KNN`,HPLC_100_out$`Missclasssified samples SVM`,HPLC_100_out$`Missclassified samples RF`),nrow = 3,byrow = T))/HPLC_100_out$`Sum of Class samples`)
colnames(HPLC_misclassified_df)=c('KNN','SVM','RF')
rownames(HPLC_misclassified_df)=c('1','2','3')

#Plot the barplots for the summary of the missclassified classification data over the 100 iterations
Enose_misclassified_plot_df=Enose_misclassified_df%>%rownames_to_column('Class')%>%pivot_longer(cols = !Class,names_to = 'Method',values_to = 'Misclass_Proportion')
ggplot(data=Enose_misclassified_plot_df)+geom_bar(mapping=aes(x=Class,y=Misclass_Proportion,fill=Method),stat = 'identity',position = 'dodge',width=.5)+ggtitle('Enose Missclassification distribution')+ylab('Missclassification Proportion')

HPLC_misclassified_plot_df=HPLC_misclassified_df%>%rownames_to_column('Class')%>%pivot_longer(cols = !Class,names_to = 'Method',values_to = 'Misclass_Proportion')
ggplot(data=HPLC_misclassified_plot_df)+geom_bar(mapping=aes(x=Class,y=Misclass_Proportion,fill=Method),stat = 'identity',position = 'dodge',width=.5)+ggtitle('HPLC Missclassification distribution')+ylab('Missclassification Proportion')


#VarImp plot for one of the iterations for each model on Enose data
Enose_varimp=grid.arrange(plot(varImp(Enose_100_out$`KNN model`),main='KNN'),plot(varImp(Enose_100_out$`SVM model`),main='SVM'),ncol=2)
Enose_varimp_RF=varImpPlot(Enose_100_out$`Random Forest Model`,main = 'RF')

#VarImp plot for one iterations for each of the model on the HPLC data
HPLC_varimp=grid.arrange(plot(varImp(HPLC_100_out$`KNN model`),main='KNN'),plot(varImp(HPLC_100_out$`SVM model`),main='SVM'),ncol=2)
HPLC_varimp_RF=varImpPlot(HPLC_100_out$`Random Forest Model`,main='RF')

#Plot Predicted v. Actual values for one of the 100 iterations of the regression models
Regression_plot_function(Enose_100_pseudo)
Regression_plot_function(Enose_100_TVC)
Regression_plot_function(HPLC_100_TVC)
Regression_plot_function(HPLC_100_pseudo)




