merge_func=function(x,y){
  return(merge(x,y,all.x=T,by='row.names')%>%column_to_rownames('Row.names'))
}

k.results=function(n,trn,tstS,trnCl,tstCl){
  acc_k=list()
  confusion_matrices=c()
  for (i in c(1:n)){
    model=knn(trn,tstS,cl = trnCl,k = i)
    confusion.matrix=confusionMatrix(model,testCl,positive = '3')
    acc_k=append(acc_k,confusion.matrix$overall[1])
    confusion_matrices=append(confusion_matrices,list(confusion.matrix))}
  return(list(acc_k,confusion_matrices))}

random_forest_tune_classification=function(data,mtry,ntree,nodesize,maxnodes,fold,nrepeats,seed){
  trcontrol=trainControl(method = 'repeatedcv',search = 'grid',number = fold,repeats = nrepeats)
  RF_grid=expand.grid(.mtry=mtry,.ntree=ntree,.nodesize=nodesize,.maxnodes=maxnodes)
  set.seed(seed)
  custom <- train(sensory~., data=data, 
                  method=customRF,metric='Accuracy', 
                  tuneGrid=RF_grid, 
                  trControl=trcontrol)
  return(custom)
}

random_forest_tune_regression=function(data,mtry,ntree,nodesize,seed){
  set.seed(seed)
  train_data=data%>%select(-ncol(.))
  response=data[,ncol(data)]
  ret_val=tune(randomForest,train.x = train_data,train.y = response,tunecontrol=ctr,ranges = list(mtry=mtry,ntree=ntree,nodesize=nodesize))
  return(ret_val$best.parameters)}

Knn_tune=function(data,k){
  train_data=data%>%select(-ncol(.))
  response=data[,ncol(data)]
  ret_val=tune(gknn,train.x = train_data,train.y=response,tunecontrol=ctr,ranges = list(k=k),scale=T)
  return(ret_val$best.parameters)}

svmrd_tune=function(data,cost,sigma,seed){
  #ctr=tune.control(sampling = 'cross',cross= fold,nrepeat= nrepeats)
  set.seed(seed)
  svm_tune=train(method='svmRadial',sensory~.,data = data,trControl=train_control,tuneGrid = expand.grid(C=cost,sigma=sigma),preProcess=c('center','scale'))
  return(svm_tune$bestTune)}

spls_tune=function(data,num_lv,ncomp,seed){
  set.seed(seed)
  predictor=data%>%select(-ncol(.))
  response=data[,ncol(data)]
  ret_val=tune.spls(X = predictor,Y=response,test.keepX = num_lv,mode = 'regression',folds = 3,nrepeat = 10,progressBar = F,ncomp = ncomp,measure='MSE',scale=T)
  return(list(ret_val$choice.keepX,ret_val$choice.ncomp$ncomp))}

confusion_matrix_extract=function(conf_matrix){
  confmatrixtable=data.frame(conf_matrix$table)
  All_sum=data.frame(confmatrixtable%>%group_by(Reference)%>%summarise(Sum=sum(Freq))%>%ungroup())
  rownames(All_sum)=All_sum$Reference
  All_sum=All_sum%>%select(-Reference)
  Incorrect_class=data.frame(confmatrixtable%>%filter(Reference!=Prediction)%>%group_by(Reference)%>%summarise(Incorrect=sum(Freq))%>%ungroup())
  rownames(Incorrect_class)=Incorrect_class$Reference
  Incorrect_class=Incorrect_class%>%select(-Reference)
  All_sum$Incorrect=Incorrect_class[rownames(All_sum),"Incorrect"]
  All_sum=All_sum%>%mutate(Proportion_incorrect=(Incorrect/Sum))
  return(All_sum)
}

extract_metrics=function(metric){
  means=c()
  std_dev=c()
  conf_int=list()
  model_names=c('KNN','RF','PLS')
  for(i in c(1:3)){
    t_test=t.test(metric[[i]])
    means=c(means,t_test$estimate[[1]])
    std_dev=c(std_dev,t_test$stderr*10)
    conf_int=append(conf_int,c(t_test$conf.int[1],t_test$conf.int[2]))
  }
  names(means)=model_names;names(std_dev)=model_names
  conf_int=data.frame(matrix(conf_int,ncol = 2,byrow = T))
  rownames(conf_int)=model_names
  ret_val=list(means,std_dev,conf_int)
  names(ret_val)=c('Mean','Std_dev','Confidence Interval')
  return(ret_val)
}


Regression_plot_function=function(output){
  par(mfrow = c(2,2))
  KNN_plot=plot(x=output$Predicted_KNN,y=output$Test_values,type='p',xlab='Predicted',ylab='Test values',main='KNN',axes=T)+abline(coef=c(0,1))+lines(x=output$Test_values+1,y=output$Test_values,col='red',lty=2)+lines(x=output$Test_values-1,y=output$Test_values,col='red',lty=2)
  RF_plot=plot(x=output$Predict_RF,y=output$Test_values,type='p',xlab='Predicted',ylab='Test values',main='RF',axes=T)+abline(coef=c(0,1))+lines(x=output$Test_values+1,y=output$Test_values,col='red',lty=2)+lines(x=output$Test_values-1,y=output$Test_values,col='red',lty=2)
  PLS_plot=plot(x=output$Predicted_PLS,y=output$Test_values,type='p',xlab='Predicted',ylab='Test values',main='PLS',axes=T)+abline(coef=c(0,1))+lines(x=output$Test_values+1,y=output$Test_values,col='red',lty=2)+lines(x=output$Test_values-1,y=output$Test_values,col='red',lty=2)
  ret_val=list(KNN_plot,RF_plot,PLS_plot)
  
  return(ret_val)
}


iter_function=function(data,k,mtry,ntree,maxnodes,nodesize,C,sigma,iter,seed){
  data=as.data.frame(data)
  knn_accuracy=c()
  RF_accuracy=c()
  svm_accuracy=c()
  class_counts=list()
  knn_miss=list()
  svm_miss=list()
  rf_miss=list()
  rand_int=sample(c(1:iter),1)
  trControl=trainControl(method = 'none')
  set.seed(seed)
  for (i in c(1:iter)){
    trainindex=createDataPartition(data$sensory,p = .7,times = 1,list=F)
    train_data=data[trainindex,]
    test_data=data[-trainindex,]
    traincl=train_data$sensory
    testcl=test_data$sensory
    train_data=train_data%>%select(-sensory)
    test_data=test_data%>%select(-sensory)
    #train_knn=preProcess(train_data,method=c('center','scale'))
    #test_knn=predict(train_knn,newdata=test_data)
    #train_knn=predict(train_knn,newdata=train_data)
    
    #Train Knn_model
    knn_fit=train(x=train_data,y=traincl,preProcess=c('center','scale'),method='knn',tuneGrid=data.frame(k=k),trControl=trControl)
    knn_predict=predict(knn_fit,test_data)
    #Train Random_forest
    random_forest_model=randomForest(x=train_data,y=traincl,ntree = ntree,maxnodes = maxnodes,mtry = mtry,nodesize = nodesize)
    random_forest_predict=predict(random_forest_model,newdata = test_data,type = 'class')
    #Train SVM
    svm.fit=train(method='svmRadial',x=train_data,y=traincl,tuneGrid=data.frame(C=C,sigma=sigma),preProcess=c('center','scale'),trControl=trControl)
    svm.predict=predict(svm.fit,test_data)
    #Confusion matrices:
    confusion_knn=confusionMatrix(knn_predict,testcl)
    confusion_rf=confusionMatrix(random_forest_predict,testcl)
    confusion_svm=confusionMatrix(svm.predict,testcl)
    #Extract inaccurate classification proportion
    knn_table=confusion_matrix_extract(confusion_knn)
    class_counts=append(class_counts,knn_table$Sum)
    knn_miss=append(knn_miss,knn_table$Incorrect)
    
    svm_table=confusion_matrix_extract(confusion_svm)
    svm_miss=append(svm_miss,svm_table$Incorrect)
    
    rf_table=confusion_matrix_extract(confusion_rf)
    rf_miss=append(rf_miss,rf_table$Incorrect)
    
    #append accuracies
    knn_accuracy=c(knn_accuracy,confusion_knn$overall[1]);RF_accuracy=c(RF_accuracy,confusion_rf$overall[1]);svm_accuracy=c(svm_accuracy,confusion_svm$overall[1])
    
    #Select Random Iteration
    if(i==rand_int){
      out_model_knn=knn_fit
      out_model_rf=random_forest_model
      out_model_svm=svm.fit
    }}
  
  #Create Matrices
  class_counts=colSums(matrix(as.numeric(class_counts),ncol=3,byrow = T))
  knn_miss=colSums(matrix(as.numeric(knn_miss),ncol = 3,byrow = T))
  svm_miss=colSums(matrix(as.numeric(svm_miss),ncol = 3,byrow = T))
  rf_miss=colSums(matrix(as.numeric(rf_miss),ncol = 3,byrow = T))
  
  ret_val=list(mean(knn_accuracy),mean(RF_accuracy),mean(svm_accuracy),class_counts,knn_miss,svm_miss,rf_miss,out_model_knn,out_model_rf,out_model_svm)
  names(ret_val)=c('Mean knn accuracy','Mean RF accuracy','Mean SVM accuracy','Sum of Class samples','Missclassified Samples KNN','Missclasssified samples SVM','Missclassified samples RF','KNN model','Random Forest Model','SVM model')
  
  return (ret_val)}


iter_function_regression=function(data,k,mtry,ntree,nodesize,ncomp,iter,seed){
  RMSE_knn=c()
  RMSE_rf=c()
  RMSE_pls=c()
  rand_int=sample(c(1:iter),1)
  set.seed(seed)
  for (i in c(1:iter)){
    trainindex=createDataPartition(data$Bac_count,times = 1,list=F)
    train_data=data[trainindex,]
    test_data=data[-trainindex,]
    train_count=train_data$Bac_count
    test_count=test_data$Bac_count
    train_data=train_data%>%select(-Bac_count)
    test_data=test_data%>%select(-Bac_count)
    train_scale=preProcess(train_data,method=c('center','scale'))
    test_scale=predict(train_scale,newdata=test_data)
    train_scale=predict(train_scale,newdata=train_data)
    
    #Train and predict KNN
    knn.fit=knnreg(x=train_scale,y=train_count,k=k)
    knn.predict=predict(knn.fit,test_scale)
    RMSE_knn=c(RMSE_knn,RMSE(knn.predict,test_count))
    
    #Train and predict Random Forest
    random_forest_model=randomForest(x=train_data,y=train_count,ntree = ntree,mtry = mtry,nodesize = nodesize)
    random_forest_predict=predict(random_forest_model,newdata = test_data)
    RMSE_rf=c(RMSE_rf,RMSE(random_forest_predict,test_count))
    
    #Train and predict pls
    pls_train=pls(X=train_data,Y = train_count,ncomp = ncomp,scale = T)
    pls_predict=predict(pls_train,newdata=test_data)
    RMSE_pls=c(RMSE_pls,RMSE(pls_predict$predict,test_count))
    
    #Extract Random Iteration
    if (i==rand_int){
      out_predict_knn=knn.predict
      out_predict_rf=random_forest_predict
      out_predict_pls=pls_predict$predict
      out_test=test_count
    }}
  
  ret_val=list(RMSE_knn,RMSE_rf,RMSE_pls,out_predict_knn,out_predict_rf,out_predict_pls,out_test)
  names(ret_val)=c('RMSE values for knn','RMSE values for Random Forest','RMSE values for PLS','Predicted_KNN','Predict_RF','Predicted_PLS','Test_values')
  return(ret_val)}




