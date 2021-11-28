rm(list=ls())
library(Ball)
library(glmnet)
library(grpreg)
expit<-function(x)
  exp(x)/(1+exp(x))


### Causal.cor
### FUNCTION:  calculate conditional ball covariance 
###            between x and y given z.
###
### INPUT:     x, a n*p matrix, samples of variable X.  
###            y, a n*q matrix, samples of variable y.
###            z, a n*1 binary vector, samples of   z . 
###            distance, if distance = TRUE, the elements
###            of x and y are considered as distance matrices.
###
### OUTPUT:    double; conditional ball covariance 
###               between x and y given z.


Causal.cor <- function(x, y, z, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha = length(index0)/length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]
    y0 <- y[index0, index0]
    x1 <- x[index1, index1]
    y1 <- y[index1, index1]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, x0, distance = TRUE)^2 + (1-alpha)*bcov(x1, y1, distance = TRUE)^2
  } else {
    x0 <- x[index0, ]
    y0 <- y[index0, ]
    x1 <- x[index1, ]
    y1 <- y[index1, ]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, y0)^2 + (1-alpha)*bcov(x1, y1)^2
  }
  return(Causal.cor)
}


### create_weights
### FUNCTION: create weight for n samples, this function is use for 
###           tuning parameter selection
###           
###
### INPUT:     ps, a n*1 vector, estimated propensity score for each sample
###            D , a n*1 vector, treatment =1 if treated; = 0 under control
###
### OUTPUT:    weight, a n*1 vector weight for each sample.


create_weights = function(ps,D){
  weight = ps^(-1)
  weight[D==0] = (1 - ps[D==0])^(-1)
  return(weight)
}


### wAMD_function
### FUNCTION: calculate weighted absolute mean difference for each choice of 
### tuning parameter. This function is use for selecting tuning parameter.
###           
###
### INPUT:     DataM   , a n*(p+2) matrix, contains both covariate, treatment and outcome.
###            varlist , a p*1  vector,names for each covariate.
###            trt.var , name of treatment, in our simulation it is "D".
###            beta    , a p*1 vector, coefficient for each covariate.
###            

### OUTPUT: a list (a) diff_vec , mean difference for each covaiate.
###                (b) wdiff_vec, weighted mean difference for each covaiate.
###                   (diff_vec*|\hat{beta}_j|)  
###                (c) wAMD     , sum of weighted mean difference
###
###
###
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

### DR double robust estimator
### FUNCTION: calculate DR estimator for each sample(which means averaged), for future analysis.
###           
###
### INPUT:     ps   a n*1  vector: estimated propensity score for each samples.
###            or1  a n*1  vector: estimated outcome regression: \hat{b}_1(X_i) for each sample.
###            or0  a n*1  vector: estimated outcome regression: \hat{b}_0(X_i) for each sample.
###            D    a n*1  vector: treatment =1 if treated; = 0 under control.
###            Y    a n*1  vector: observed outcome.
###
### Output:    n*1 DR estimates   
###
###

DR =  function(ps,or1,or0,D,Y){
  (D*Y-(D-ps)*or1)/ps -  (   (1-D)*Y +(D-ps)*or0)/(1-ps)
}


### CBS
### FUNCTION: Proposed CBS method for estimating average treatment effect.
###           
###
### INPUT:     
###            X:   a n*p  matrix: n*p covariate.
###            D    a n*1  vector: treatment =1 if treated; = 0 under control.
###            Y    a n*1  vector: observed outcome.
###            alpha: double; significant level for confidence interval construction.
      
### OUTPUT: a list (a) point estimate, the doubly debiased lasso estimator.
###                (b) lower bound, lower bound of the interval.
###                (b) upper bound, upper bound of the interval.
###                (c) Variance, a double; the variance of the dblasso estimator.
###
###
###

CBS <- function(X,D,Y,alpha = 0.05){
  p = ncol(X)
  n = nrow(X)
  threshold = min(30,p) 
  #set threshold for screening.
  ballcor<-rep(NA, p)
  for (j in 1:p){
    # calculate conditional ball covariance for each variable.
    ballcor[j]<-Causal.cor(X[,j],Y,D)
  }
  
  
  var.list_Ball = c(paste("X",1:threshold,sep=""))
  var.list = c(paste("X",1:p,sep=""))
  
  # set lambda values for grid search.
  
  lambda = log(max(p,n))^0.75/n^0.5
  lambda_val = seq(lambda*0.1,lambda*10,length.out= 10)
  names(lambda_val) = paste0("lambda",1:10)
  
  # set gamma values for grid search.
  
  gamma_vals = seq(3,20,1)
  names(gamma_vals) = paste0("gamma",gamma_vals)
  
  
  # screening procedure 
  ballorder<-order(ballcor, decreasing=TRUE)

  # select the top threshold number of variables
  ballsisselectedindex<-ballorder[1:threshold]
  ballsisselectedindex = ballsisselectedindex[order(ballsisselectedindex)]
  weight = ballcor[ballsisselectedindex]
  
  # the data matrix after screening
  Data = X[,ballsisselectedindex]
  Data = as.data.frame(Data)
  names(Data) = var.list_Ball
  Data$D = D
  Data$Y = Y
  
  
  # centerlize, standardlize
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # weight for each variable for refined selection
  betaXY = weight
  betaXY = weight/max(weight)
  
  ## Want to save wAMD and propensity score coefficients for
  ## each lambda and gamma value
  
  wAMD_vec = rep(NA, length(lambda_val)*length(gamma_vals))
  names(wAMD_vec) = paste0( rep(names(lambda_val),each = length(gamma_vals)),names(gamma_vals))
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(wAMD_vec)))
  names(coeff_XA) = names(wAMD_vec)
  rownames(coeff_XA) = var.list_Ball
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda and gamma value ####################
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  
  
  for( lil in names(lambda_val) )
    for(mim in names(gamma_vals)){
      il = lambda_val[lil]
      ig = gamma_vals[mim]
      fitx = as.matrix(Data[,1:threshold])
      alasso <- glmnet(x = fitx, y = Data$D,
                        type.measure = "class",
                        family = "binomial",
                        alpha = 1,
                        penalty.factor = c(abs(betaXY)^(-ig)),
                        lambda = il)
    
    # calculate propensity score 
    prob = predict(alasso,newx = fitx)
    prob = exp(prob)/(1+exp(prob))
    Data[,paste("f.pA",paste0( lil,mim),sep="")] = prob
    # save propensity score coefficients
    coeff_XA[var.list_Ball,paste0( lil,mim)] = coef.glmnet(alasso,s = alasso$lambda.min)[var.list_Ball,]
    # create inverse probability of treatment weights
    Data[,paste("w",paste0( lil,mim),sep="")] = create_weights(ps=Data[,paste("f.pA",paste0( lil,mim),sep="")],D=Data$D)
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[paste0( lil,mim)] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="D",
                                  wgt=paste("w",paste0( lil,mim),sep=""),beta=betaXY)$wAMD
    # save ATE estimate for this lambda value
  } # close loop through lambda values
  
  # find the target (gamma,lambda) with smallest wAMD score.
  tt = which.min(wAMD_vec)
  
  #save the estimated propensity score model
  fitted.ps = Data[,paste("f.pA",names(tt),sep="")]
  #outcome regression
  {
    # use lasso to fit the treated group.
    X1 = X[D==1,]
    Y1 = Y[D==1]
    X0 = X[D==0,]
    Y0 = Y[D==0]
    
    lasso<-cv.glmnet(x = X1, y = Y1,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso1 <- coef.glmnet(lasso,lasso$lambda.min)
    index1 = which(coeflasso1[-1]!=0)
    X1.fit = X1[,index1]
    
    data1 = data.frame(Y1,X1.fit)
    names(data1) = c("Y",paste0("v",1:length(index1)))[1:ncol(data1)]
    
    ormodel1 <- lm(Y~.,data1)
    
    data.fit = data.frame(Y,X[,index1])
    names(data.fit) = c("Y",paste0("v",1:length(index1))) [1:ncol(data1)]
    # save predict for the treated group.
    orfit1 = predict.lm(ormodel1, newdata = data.fit)
    
    # use lasso to fit the controlled group.
    lasso<-cv.glmnet(x = X0, y = Y0,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso0 <- coef.glmnet(lasso,lasso$lambda.min)
    index0 = which(coeflasso0[-1]!=0)
    X0.fit = X0[,index0]
    
    data0 = data.frame(Y0,X0.fit)
    names(data0) = c("Y",paste0("v",1:length(index0)))[1:ncol(data0)]
    
    ormodel0 <- lm(Y~.,data0)
    data.fit = data.frame(Y,X[,index0])
    names(data.fit) = c("Y",paste0("v",1:length(index0))) [1:ncol(data0)]
    #save the estimated data for the controlled group. 
    orfit0 = predict.lm(ormodel0, newdata = data.fit)
  }
  # get the double robust estimate.
  result = DR(fitted.ps,orfit1,orfit0,D,Y) 
  
  # get the point estimate and variance of the resulting estimator. 
  result = c(mean(result),var(result))
  c("point estimate" = result[1],
    "lower bound" = result[1]-qnorm(1-alpha/2)*sqrt(result[2]/n),
    "upper bound" = result[1]+qnorm(1-alpha/2)*sqrt(result[2]/n),
    "variance" =  sqrt(result[2]/n))
}

#generate data
set.seed(12345*2021+54321)
n = 600
p = 200
X = matrix(runif(n*p,min = -1, max = 1),ncol = p)
prob = rowSums(X[,1:2])*0.2 + rowSums(X[,5:6])*0.3
prob = exp(prob)/(1+exp(prob))
D = rbinom(n,size = 1 , prob)
Y = D*2 + rowSums(X[,c(1:4)])*2 + rnorm(n)

CBS(X,D,Y)
  