library(dismo)
library(gbm)
library(pROC)
library(tidyverse)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
library(sf)
library(ggspatial)
library(PRROC)
library(factoextra)
library(readxl)

setwd('')
var<-read.csv('var.csv') # independent variables
var_name<-read.csv('varname.csv')[,1] # variable names
case_all <- read.csv('rat_fenbu_2.csv')# rat data used to model
moi_study<-names(case_all)[-1]

###logistic model to calculate the weight of sampleing for reducing investigation error----
case_all$if_survey <- ifelse(case_all$ADCODE99 %in% all_survey,1,0)
case_logistic<- var[,c('ADCODE99',var_name)] %>%   
  dplyr::left_join(case_all[,c('ADCODE99','if_survey')],by=c('ADCODE99'))
aa <- glm(as.formula(str_c('if_survey~',str_c(var_name,collapse = '+'))),
          data=case_logistic,family = 'binomial')
bb <- step(aa,direction = 'backward')
dd <- glm(as.formula(str_c('if_survey~',
                           str_c(names(bb$coefficients)[-which(names(bb$coefficients)%in%c("(Intercept)",'value_43','value_45','value_46','gdp'))],collapse = '+'))),
          data=case_logistic,family = 'binomial')
a <- predict(dd,newdata=case_logistic,type='response')
case_var$weight <- 1/a
case_var$weight <- case_var$weight/mean(case_var$weight[case_logistic$if_survey==1])

### built model
for(m_big_i in 1:46){
  ###m=targeted rat, case_var_model= data including all variables
  setwd('') #each root to spefic mosquito
  m<-moi_study[m_big_i]
  id=which(names(case_var)=='ADCODE99')
  var_c=which(names(case_var) %in% var_name)
  var_mo=which(names(case_var)== m)
  ir_inf<-10;ir_mo<-100
  
  data<-case_var[,c(var_mo,id,var_c,which(names(case_var)=='weight'))]
  names(data)[1]<-'case'
  var_model<-which(names(data) %in% var_name)
  
  data1<-data[data$case==1&!is.na(data$case),]
  data0<-data[data$case==0&!is.na(data$case),]
  n1<-ceiling(0.75*nrow(data1))
  n0<-ceiling(0.75*nrow(data0))
  dir.create(m)
  setwd(m)
  write.csv(rbind(data1,data0),'all.csv',row.names = F)
  write.csv(data,'all_NA.csv',row.names = F)
  for (i in 1:150){ #randomly split data for 150 times
    set.seed(1234+i)
    train1<-sample(1:nrow(data1),n1)
    set.seed(1234+i)
    train0<-sample(1:nrow(data0),n0)
    dftrain1<-data1[train1,];dftest1<-data1[-train1,];dftrain0<-data0[train0,];dftest0<-data0[-train0,]
    dftrain<-rbind(dftrain1,dftrain0);dftest<-rbind(dftest1,dftest0)
    dir.create('train');traini<-paste0('train\\train',i,'.csv')
    dir.create('test');testi<-paste0('test\\test',i,'.csv')
    write.csv(dftrain,traini,row.names = F);write.csv(dftest,testi,row.names = F)
  }
  ###step2.run models for 10 times with all variables to calculate relative influence
  model.all<-list();inf<-data.frame(var=c(names(data)[var_model]))
  while (ncol(inf)<(ir_inf+1)){
    i<-ncol(inf)
    train.data<-read.csv(paste0('train\\train',i,'.csv'))
    modeli<-gbm.step(data=train.data,gbm.x = var_model,gbm.y = 1,family = 'bernoulli', tree.complexity = 5,
                     site.weights =train.data$weight,
                     learning.rate = 0.005, bag.fraction = 0.75,n.folds = 10,max.trees = 3000,silent = T, plot.main = F)
    model.all[[i]]<-modeli
    if(is.null(modeli)) next()
    infi<-summary(modeli);names(infi)<-c('var',paste0('inf',i))
    inf<-left_join(inf,infi,by='var')
    print(paste0(m,i))
  }
  inf <- inf %>%
    mutate(mean=apply(inf[,str_c('inf',1:10)],1,mean))
  
  var_inf<-function(inf,cut=2){
    aa<-as.character(inf[inf$mean>cut,1])
    a<-c(seq_along(data))[names(data) %in% aa]
    a
  }
  var_m<-var_inf(inf,cut = 2)
  
  ###step.3 run models for 100 times with filtered variables as final models
  model_z<-list();ROC<-list()
  i<-1;n_mo<-10;id_model_z<-1
  c_train<-c();c_test<-c()
  while(length(model_z)<(ir_mo+1)){
    n_traini<-paste0('train\\train',i+n_mo,'.csv');n_testi<-paste0('test\\test',i+n_mo,'.csv');
    train.data<-read.csv(n_traini);test.data<-read.csv(n_testi);all.data<-bind_rows(train.data,test.data)
    set.seed(12345+i)
    model_zi<-gbm.step(data<-train.data,gbm.x = var_m,gbm.y = 1,family = 'bernoulli', tree.complexity = 5,
                       site.weights =train.data$weight,
                       learning.rate = 0.005, bag.fraction = 0.75,n.folds = 10,max.trees = 3000,silent = T, plot.main = F)
    print(str_c(i,'_100'))
    if(!is.null(model_zi)) {model_z[[id_model_z]]<-model_zi;id_model_z<-id_model_z+1
    c_train<-c(c_train,n_traini);c_test<-c(c_test,n_testi)}
    i<-i+1
  }
  pre_df<- function(name_df){
    df_re<-NULL
    for(i in 1:100){
      if(length(name_df)==100) data<-read.csv(name_df[i]) else data<-read.csv(name_df)
      codei<-str_c('code',i);orii<-str_c('ori',i);predi<- str_c('pre',i)
      preds_i<- predict.gbm(model_z[[i]],data,n.trees = model_z[[i]]$n.trees,type='response')
      if(is.null(df_re)) df_re<-cbind(data[,c(2,1)],preds_i) else df_re<-cbind(df_re,data[,c(2,1)],preds_i)
      names(df_re)[(i*3-2):(i*3)]<-c(codei,orii,predi)
    }
    df_re   
  }
  df_pred_train<-pre_df(c_train[1:100])
  df_pred_test<-pre_df(c_test[1:100])
  df_pred_all<-pre_df('all.csv')
  df_pred_allna<-pre_df('all_NA.csv')
  
  
  roc_wt<-function(df_pred){
    roc_l<-list();auc<-c();cut<-c()
    for(i in 1:100){
      roc_l[[i]]<-roc(df_pred[,(i*3-1)],df_pred[,(i*3)])
      auc<-c(auc,roc_l[[i]]$auc)
      cut<-c(cut,coords(roc_l[[i]],'best',transpose=T)['threshold'])
    }
    roc_l[[101]]<- auc;roc_l[[102]]<- cut
    roc_l
  }
  roc_train<-roc_wt(df_pred_train);roc_test<-roc_wt(df_pred_test);roc_all<-roc_wt(df_pred_all)
  
  
  df_roc<-data.frame(name=c('train','test','all'),
                     auc=c(mean(roc_train[[101]]),mean(roc_test[[101]]),mean(roc_all[[101]])),
                     cut=c(mean(roc_train[[102]]),mean(roc_test[[102]]),mean(roc_all[[102]])))
  
  
  df_pred_allna$mean<-apply(df_pred_allna[,(1:100)*3],1,mean)
  df_pred_allna$if_not<-ifelse(df_pred_allna$mean>mean(roc_all[[102]],na.rm = T),1,0)
  if_ci<-function(a){
    a_n<-c()
    for(i in 1:100){
      if(is.na(roc_all[[102]][i])) ai<-NA else {if(a[i*3]>roc_all[[102]][i]) ai<-1 else ai<-0}
      a_n<-c(a_n,ai)
    }
    sum(a_n,na.rm = T)/sum(!is.na(a_n))
  }
  df_pred_allna$if_n<- apply(df_pred_allna,1,if_ci)
  df_pred_allna$if_n_pan<- ifelse(df_pred_allna$if_n>0.5,1,0)
  save.image()
  ####save result
  dir.create('result')
  write.csv(df_pred_train,'result\\pred_train.csv')
  write.csv(df_pred_test,'result\\pred_test.csv')
  write.csv(df_pred_all,'result\\pred_all.csv')
  write.csv(df_pred_allna,'result\\pred_allna.csv')
  write.csv(df_roc,'result\\roc.csv')
}


###Part 2. GBRT for diseases
library(xgboost)
library(readxl)
library(tidyverse)
library(rBayesianOptimization)
library(ROCR)
setwd('')
fading<-read.csv('incidence.csv') #Statutory reported incidence data
pop<-read_xlsx('code.xlsx')[,c('ADCODE99','pop_2010')] #population
collect1<-read_xlsx('collect.xlsx') #data from mosquito and host animal
v_infor <- read_xlsx('disease.xlsx') #Name of disease
class <- read_xlsx('class.xlsx') #Class of disease
#Data processing
data_list <- list()
for(d in 1:2){
  aa_f <- fading %>% filter(disease==v_infor$disease[d]) %>% 
    group_by(code) %>% 
    summarise(n=sum(n))
  aa_c <- collect1 %>% filter(disease %in% class$species[class[v_infor$col[d]]==v_infor$coll[d]])
  aa_f$code <- as.numeric(aa_f$code)
  aa_d <- pop
  aa_d <- aa_d %>% left_join(aa_f,by=c('ADCODE99'='ADCODE99')) 
  aa_d$rat <- ifelse(aa_d$ADCODE99 %in% aa_c$ADCODE99,1,0)
  aa_d$exist <- ifelse(aa_d$ADCODE99 %in% aa_c$ADCODE99|!is.na(aa_d$n),1,0)
  aa_d$incidence <- aa_d$n/aa_d$pop_2010/15*100000
  data_list[[d]] <- aa_d
}
names(data_list) <- c('Hantaviridae','Leptospira')


##var_的蚊子合并----
var <- read.csv('var_0118.csv') #The independent variables

##built model
name_all<-c('Hantaviridae','Leptospira')

for(id in 1:2){
  name<-name_all[id]
  setwd('')
  dir.create(name)
  var_name <- intersect(var_list$Var[!is.na(var_list[name])],names(var))
  data_id <- data_list[[id]][c('ADCODE99','exist','incidence')] %>% left_join(var[c('ADCODE99',var_name)],by=c('ADCODE99'))
  ##part.1 logistic regression
  params<-list(eta=0.005,max_depth=8,subsample=0.75,
               objective = "reg:logistic",eval_metric = "auc")
  data_id1<-data_id %>%
    filter(exist==1)
  data_id0<-data_id%>%
    filter(exist==0)
  inf_t<-NULL
  ##run models for 10 times with all variables to calculate relative influence
  for(i in 1:10){
    set.seed(i+1234);tt1<-sample(1:nrow(data_id1),round(0.75*nrow(data_id1)))
    set.seed(i+1234);tt0<-sample(1:nrow(data_id0),round(0.75*nrow(data_id0)))
    dtrain<-rbind(data_id1[tt1,],data_id0[tt0,]);dtest<-rbind(data_id1[-tt1,],data_id0[-tt0,])
    x_train<-xgb.DMatrix(as.matrix(dtrain[-c(1:3)]),label=as.matrix(dtrain[2]))
    x_test<-xgb.DMatrix(as.matrix(dtest[-c(1:3)]),label=as.matrix(dtest[2]))
    watchlist<-list(train=x_train,eval=x_test)
    set.seed(i+1234)
    bst_cv <- xgb.cv(params = params,data=x_train,watchlist = watchlist,verbose = 1,nfold = 5,nround = 100,
                     early_stopping_rounds = 10)
    bst<-xgb.train(params = params,data=x_train,nrounds = bst_cv$best_iteration,watchlist = watchlist,verbose = 1)
    inf<-xgb.importance(model=bst)
    if(is.null(inf_t)) inf_t<-inf[,1:2] else inf_t<-full_join(inf_t,inf[,1:2],by='Feature')
  }
  inf_t$Gain<-apply(inf_t[,2:11],1,mean,na.rm=T)
  var_inf<-inf_t$Feature[inf_t$Gain>=0.015&!is.na(inf_t$Gain)]
  #run models for 100 times with filtered variables 
  df_yn_f<-data_id[c('exist',var_inf)]
  x_pred<-xgb.DMatrix(as.matrix(df_yn_f[-1]),
                      label=as.matrix(df_yn_f[1]))
  df_yn1<-df_yn_f %>%
    filter(exist==1)
  df_yn0<-df_yn_f%>%
    filter(exist==0)
  
  inf_f<-NULL
  auc_train<-c();auc_test<-c()
  bst_list<-list();pred_f<-df_yn_f['exist']
  t_all <- c();b_all <- 0
  for(i in 1:100){
    set.seed(i+1234);tt1<-sample(1:nrow(df_yn1),round(0.75*nrow(df_yn1)))
    set.seed(i+1234);tt0<-sample(1:nrow(df_yn0),round(0.75*nrow(df_yn0)))
    dtrain<-rbind(df_yn1[tt1,],df_yn0[tt0,]);dtest<-rbind(df_yn1[-tt1,],df_yn0[-tt0,])
    x_train<-xgb.DMatrix(as.matrix(dtrain[-1]),label=as.matrix(dtrain[1]))
    x_test<-xgb.DMatrix(as.matrix(dtest[-1]),label=as.matrix(dtest[1]))
    watchlist<-list(train=x_train,eval=x_test)
    set.seed(i+1234)
    bst_cv <- xgb.cv(params = params,data=x_train,watchlist = watchlist,verbose = 1,nfold = 5,nround = 100,
                     early_stopping_rounds = 10)
    
    bst<-xgb.train(params = params,data=x_train,nrounds = bst_cv$best_iteration,watchlist = watchlist,verbose = 1)
    bst_list[[i]]<-bst
    inf<-xgb.importance(model=bst)
    if(is.null(inf_f)) inf_f<-inf[,1:2] else inf_f<-full_join(inf_f,inf[,1:2],by='Feature')
    auc_train<-c(auc_train,as.numeric(bst$evaluation_log[nrow(bst$evaluation_log),2]))
    auc_test<-c(auc_test,as.numeric(bst$evaluation_log[nrow(bst$evaluation_log),3]))
    pred<-predict(bst,newdata = x_pred)
    pred_f<-cbind(pred_f,pred)
    df_roc<-prediction(pred,as.integer(pred_f$exist))
    tt<-performance(df_roc,measure="sens", x.measure="spec")
    t<-mean(tt@alpha.values[[1]][(tt@y.values[[1]]+tt@x.values[[1]])==max(tt@y.values[[1]]+tt@x.values[[1]],na.rm = T)],na.rm=T)
    b <- as.numeric(pred>t)
    b_all <- b_all+b
    t_all <- c(t_all,t)
  }
  pred_f$mean<-apply(pred_f[,2:101],1,mean,na.rm=T)
  inf_f$mean<-apply(inf_f[,2:101],1,mean,na.rm=T)
  qq<-pROC::roc(as.integer(pred_f$exist),pred_f$mean)
  t_mean <- mean(t_all)
  pred_f$presense<-as.integer(pred_f$mean>t_mean)

  ####Part.2 gamma regression
  params<-list(eta=0.005,max_depth=8,subsample=0.75,
               objective = "reg:gamma",
               tree_method='hist',max_bin =1500)
  df_gamma_yn<-data_id[!is.na(data_id$incidence),]
  df_gamma_yn$incidence <- df_gamma_yn$incidence*bbb[id]
  df_gamma_yn <- df_gamma_yn %>% left_join(code[c('ADCODE99','yy')],by='ADCODE99')
  df_gamma_yn$weight <- ifelse(df_gamma_yn$yy>33,0.25,1)
  df_gamma_yn <- df_gamma_yn[,c(1:3,ncol(df_gamma_yn),4:(ncol(df_gamma_yn)-2))]
  inf_t2<-NULL
  #run models for 10 times with all variables to calculate relative influence
  for(i in 1:10){
    set.seed(i+123);tt<-sample(1:nrow(df_gamma_yn),round(0.75*nrow(df_gamma_yn)))
    dtrain<-df_gamma_yn[tt,];dtest<-df_gamma_yn[-tt,]
    x_train<-xgb.DMatrix(as.matrix(dtrain[-c(1:4)]),label=as.matrix(dtrain[3]))
    x_test<-xgb.DMatrix(as.matrix(dtest[-c(1:4)]),label=as.matrix(dtest[3]))
    watchlist<-list(train=x_train,eval=x_test)
    set.seed(i+123)
    bst_cv <- xgb.cv(params = params,data=x_train,watchlist = watchlist,verbose = 1,nfold = 5,nround = 100,
                     early_stopping_rounds = 10,weight=dtrain$weight)
    set.seed(i+123)
    bst<-xgb.train(params = params,data=x_train,nrounds = bst_cv$best_iteration,watchlist = watchlist,verbose = 1,weight=dtrain$weight)
    inf<-xgb.importance(model=bst)
    if(is.null(inf_t2)) inf_t2<-inf[,1:2] else inf_t2<-full_join(inf_t2,inf[,1:2],by='Feature') 
  }
  inf_t2$mean<-apply(inf_t2[,c(2:11)],1,mean,na.rm=T)
  var_inf<-inf_t2$Feature[inf_t2$mean>=0.015]
  
  inf_t2[,c('Feature','mean')]
  ####run models for 100 times with filtered variables 
  df_yn_f2<-df_gamma_yn[,c('incidence','weight',var_inf)]
  
  data_id$incidence <- ifelse(is.na(data_id$incidence),0,data_id$incidence)
  x_pred<-xgb.DMatrix(as.matrix(data_id[var_inf]),
                      label=as.matrix(data_id['incidence']))
  inf_f2<-NULL
  bst_list2<-list();pred_f2<-data_id['incidence']
  loss_train<-c();loss_test<-c()
  for(i in 1:100){
    set.seed(i+1234);tt<-sample.int(n=nrow(df_yn_f2),size=round(0.75*nrow(df_yn_f2)),prob = df_yn_f2$weight)
    dtrain<-df_yn_f2[tt,];dtest<-df_yn_f2[-tt,]
    x_train<-xgb.DMatrix(as.matrix(dtrain[-(1:2)]),label=as.matrix(dtrain[1]))
    x_test<-xgb.DMatrix(as.matrix(dtest[-(1:2)]),label=as.matrix(dtest[1]))
    watchlist<-list(train=x_train,eval=x_test)
    set.seed(i+1234)
    bst_cv <- xgb.cv(params = params,data=x_train,watchlist = watchlist,verbose = 1,nfold = 5,nround = 100,
                     early_stopping_rounds = 10,weight=dtrain$weight)
    set.seed(i+1234)
    bst<-xgb.train(params = params,data=x_train,nrounds =bst_cv$best_iteration,watchlist = watchlist,weight=dtrain$weight)
    bst_list2[[i]]<-bst
    inf<-xgb.importance(model=bst)
    if(is.null(inf_f2)) inf_f2<-inf[,1:2] else inf_f2<-full_join(inf_f2,inf[,1:2],by='Feature') 
    pred<-predict(bst,newdata = x_pred)
    pred_f2<-cbind(pred_f2,pred)
    loss_train<-c(loss_train,as.numeric(bst$evaluation_log[20,2]))
    loss_test<-c(loss_test,as.numeric(bst$evaluation_log[20,3]))
  }
  inf_f2$mean<-apply(inf_f2[,c(2:101)],1,mean,na.rm=T)
  pred_f2$mean<-apply(pred_f2[,2:101],1,mean) 
  
  data_id[,'pre1']<-pred_f$mean
  data_id[,'pre2']<-pred_f2$mean 
  
  data_id$pred_12<-ifelse(data_id$pre1<=t_mean,0,
                          data_id$pre1*data_id$pre2)
  # Save results
  write.csv(pred_f2,paste0('鼠疾病预测/',name,'\\pred_gamma.csv'),row.names=F)
  write.csv(inf_f2,paste0('鼠疾病预测/',name,'\\inf_gamma.csv'),row.names=F)
  write.csv(pred_f,paste0('鼠疾病预测/',name,'\\pred_logistic.csv'),row.names=F)
  write.csv(inf_f,paste0('鼠疾病预测/',name,'\\inf_logistic.csv'),row.names=F)
  write.csv(data_id,paste0('鼠疾病预测/',name,'\\pred_incidence.csv'),row.names=F)
  save.image(file = paste0('鼠疾病预测/',name,'\\xgb.Rdata'))
}

