library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(caret)

rf_nodesize <- 5
seed <- 1234


#---------------- random survival forest (RSF) ----------------#


set.seed(seed)
model <- rfsrc(Surv(BCR.time,BCR)~.,data = na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]),
             ntree = 80,nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(model$err.rate)

score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=predict(model,newdata = x)$predicted)})
index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- 'RSF'


#---------------- elastic network (Enet) ----------------#


for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  model = cv.glmnet(as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,colnames(training_cohort)[-c(1:3)]]), as.matrix(Surv(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time,na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR)),family = "cox",alpha=alpha,nfolds = 10)
  score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=as.numeric(predict(model,type='link',newx=as.matrix(x[,-c(1,2)]),s=model$lambda.min)))})
  
  index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  index$Model <- paste0('Enet','[α=',alpha,']')
}


#---------------- Lasso ----------------#


set.seed(seed)
model = cv.glmnet(as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,colnames(training_cohort)[-c(1:3)]]), as.matrix(Surv(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time,na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR)),family = "cox",alpha=1,nfolds = 10)
coef.min = coef(model, s = "lambda.min") 
score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}), function(x) {
  x[, 1:2] %>% 
    mutate(RS = apply(as.matrix(x[,-c(1,2)]), 1, function(x_new) sum(coef.min[-1] * x_new)))
})
index <- data.frame(Cindex = sapply(score, function(x) as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1]))) %>% 
  rownames_to_column('ID')
index$Model <- paste0('Lasso')


#---------------- stepwise Cox ----------------#


for (direction in c("both", "backward", "forward")) {
  model <- step(coxph(Surv(BCR.time,BCR)~.,na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])),direction = direction)
  score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=predict(model,type = 'risk',newdata = x))})
  
  index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  index$Model <- paste0('StepCox','[',direction,']')
}


#---------------- CoxBoost ----------------#


na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]) <- as.data.frame(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]))
set.seed(seed)
pen <- optimCoxBoostPenalty(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR.time'],na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR'],as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR.time'],na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR'],as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
model <- CoxBoost(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR.time'],na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,'BCR'],as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=as.numeric(predict(model,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- paste0('CoxBoost')


#---------------- partial least squares regression for Cox (plsRcox) ----------------#


set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,colnames(training_cohort)[-c(1:3)]],time=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time,status=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR),nt=10,nfold = 10,verbose = F)
model <- plsRcox(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,colnames(training_cohort)[-c(1:3)]],time=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time,event=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR,nt=as.numeric(cv.plsRcox.res[5]))
score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=as.numeric(predict(model,type="lp",newdata=x[,-c(1,2)])))})

index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- paste0('plsRcox')


#---------------- supervised principal components (SuperPC) ----------------#


data <- list(x=t(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[,-c(1,2)]),y=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time,censoring.status=na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR,featurenames=colnames(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]))[-c(1,2)])
set.seed(seed)
model <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.model <- superpc.cv(model,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2, 
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$BCR.time,censoring.status=w$BCR,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(model,data,test,threshold = cv.model$thresholds[which.max(cv.model[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- paste0('SuperPC')


#---------------- generalized boosted regression modeling (GBM) ----------------#


set.seed(seed)
model <- gbm(formula = Surv(BCR.time,BCR)~.,data = na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]),distribution = 'coxph',
           n.trees = 10,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.1,
           cv.folds = 10,n.cores = 6)
best <- which.min(model$cv.error)
set.seed(seed)
model <- gbm(formula = Surv(BCR.time,BCR)~.,data = na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]),distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.1,
           cv.folds = 10,n.cores = 8)
score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=as.numeric(predict(model,x,n.trees = best,type = 'link')))})

index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- paste0('GBM')


#---------------- survival support vector machine (survival-SVM) ----------------#


model = survivalsvm(Surv(BCR.time,BCR)~., data= na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]), gamma.mu = 1)

score <- lapply(lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}),function(x){cbind(x[,1:2],RS=as.numeric(predict(model, x)$predicted))})

index <- data.frame(Cindex=sapply(score,function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
index$Model <- paste0('survival-SVM')


#---------------- Ridge ----------------#


set.seed(seed)
model <- cv.glmnet(as.matrix(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])[, colnames(training_cohort)[-c(1:3)]]),as.matrix(Surv(na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR.time, na.omit(training_cohort[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])])$BCR)), family = "cox", alpha = 0, nfolds = 10)
score <- lapply (lapply(all_cohort,function(x){x[,c('BCR.time','BCR',colnames(training_cohort)[-c(1:3)])]}), function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(model, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = model$lambda.min)))
})
index <- data.frame(Cindex = sapply(score, function(x) {
  as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])
})) %>%
  rownames_to_column('ID')
index$Model <- 'Ridge'

