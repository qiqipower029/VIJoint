rm(list=ls())
require(survival)
require(splines)
require(survival)
require(JM)
require(Matrix)
require(statmod)
require(MASS)
require(openxlsx)

#-------------------------Scenario 1-------------------------#
# n = 500, 1500, 5000
# K = 1
# Sigma = 0.5^2I


sig = 0.1 # sd for epsilon
Ngene = 1 # number of genes
Ninfo = 1 # signal=Ngene, no noise
eff = 1  # gene coefficients, true value of gam.y=1, gam.v=1
n = 500 # number of subjects
lam = 3 # Weibull parameter
theta = 10 # Weibull parameter
Sig = 0.5 # sd for MVN Sigma

Sigma = Sig^2*diag(2) # variance-covariance matrix of betas for each subject
k = 1
while (k<=100) {

AB1.list = lapply(1:n,function(i){
  matrix(mvrnorm(n = 1, rep(1.5, 2), Sigma = Sigma), nrow = Ngene)
})


Y = rep(0,n)
x = runif(n,-1,1)
beta = 1

tryCatch({
for(i in 1:n){
  AB1 = AB1.list[[i]][1,]
  Atrue = AB1[[1]]
  Btrue = AB1[[2]]
  
  logr = log(runif(1))
  ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])} # hazard function h(t)
  f = function(x){ integrate(ff,0,x)$value + logr } # Hazard function H(t)
  Y[i] = uniroot(f,c(0,50))$root
}


C = rweibull(n,2,3)
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)


len = 0.2 # the interval between each time points, also known as delta
#S = sapply(1:n,function(i){c(0,sort(runif(mmm-1,0,1))) })
data = lapply(1:n,function(i){
  if(data.id$ftime[i] <= len) {
  ss = c(0, len)
  ss = jitter(ss)
  ss[1]=0
  } else {
    ss=seq(0,data.id$ftime[i],by=len)
    ss = jitter(ss)
    ss[1]=0
}
  years = rep(ss,Ngene)
  item = rep(paste("gene",1:Ngene,sep=""),each=length(ss))
  ABtmp = AB1.list[[i]][rep(1:Ngene,each=length(ss)),,drop = F]
  
  data.frame(ID=i,item=item,years=years,AB=ABtmp)
})

data = do.call(rbind,data)
value = data[,'AB.1'] + data[,'AB.2']*data[,'years'] + rnorm(nrow(data),0,sig)
data = cbind(data,value)

SS = data.id[match(data$ID,data.id$ID),'ftime']
data = data[data$years<SS,]
}, error = function(e){})

LongData = data
SurvData = data.id

write.xlsx(LongData, file = paste0('longdata_s1',k,'.xlsx'))
write.xlsx(SurvData, file = paste0('survdata_s1',k,'.xlsx'))

# JM and mjoint data generation
marker.name = sort(unique(LongData$item))
WideData=list()
for (i in 1:n){
  tmp=LongData[LongData$ID==i,]
  nrm=length(unique(tmp$years))
  # long data
  wideData=matrix(0,nrm,length(marker.name))
  rest=tmp[tmp$item==marker.name[1],c(1,3:5)] # ID, years,AB.1, AB.2
  # surv data
  wideSurv=SurvData[SurvData$ID==i,]
  wideSurvData=matrix(rep(as.matrix(wideSurv),each=nrm),nrow=nrm)
  for (j in 1:length(marker.name)){
    wideData[,j]=tmp[tmp$item==marker.name[j],"value"]
  }
  colnames(wideData)=marker.name
  colnames(wideSurvData)=colnames(SurvData) # ID, fstat, ftime, x
  WideData[[i]]=cbind(rest,wideData,wideSurvData)
}
data.mjoint = data.frame(do.call(rbind,WideData))
save(data.mjoint, file = paste0('mjoint_s1',k, '.RData'))
k = k+1
}



#-------------------------Scenario 2-------------------------#
# n = 500
# K = 3, 5, 8
# Sigma = 0.5^2I

sig = 0.1 # sd for epsilon
Ngene = 8 # number of genes
Ninfo = Ngene # signal=Ngene, no noise
eff = 1  # gene coefficients, true value of gam.y=1, gam.v=1
n = 500 # number of subjects
lam = 3 # Weibull parameter
theta = 10 # Weibull parameter
Sig = 0.5 # sd for MVN Sigma

Sigma = Sig^2*diag(2) # variance-covariance matrix of betas for each subject
k = 1
while (k<=100) {

AB1.list = lapply(1:n,function(i){
  matrix(mvrnorm(n = Ngene, rep(1.5, 2), Sigma = Sigma), nrow = Ngene)
}) # generate random intercept and slope for each gene and each subjects


Y = rep(0,n)
x = runif(n,-1,1)
beta = 1

tryCatch({
for(i in 1:n){
  AB1 = AB1.list[[1]][1:Ninfo,]
  Atrue = sum(AB1[,1]) # summation of random intercept
  Btrue = sum(AB1[,2]) # summation of random slope
  
  logr = log(runif(1))
  ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])} # hazard function h(t)
  f = function(x){ integrate(ff,0,x)$value + logr } # Hazard function H(t)
  Y[i] = uniroot(f,c(0,50))$root
}

C = rweibull(n,2,3)
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)


len = 0.2 # the interval between each time points, also known as delta
#S = sapply(1:n,function(i){c(0,sort(runif(mmm-1,0,1))) })
data = lapply(1:n,function(i){
  
  ss=seq(0,data.id$ftime[i],by=len)
  ss = jitter(ss)
  ss[1]=0
  
  years = rep(ss,Ngene)
  item = rep(paste("gene",1:Ngene,sep=""),each=length(ss))
  ABtmp = AB1.list[[i]][rep(1:Ngene,each=length(ss)),]
  
  data.frame(ID=i,item=item,years=years,AB=ABtmp)
})

data = do.call(rbind,data)
value = data[,'AB.1'] + data[,'AB.2']*data[,'years'] + rnorm(nrow(data),0,sig)
data = cbind(data,value)

SS = data.id[match(data$ID,data.id$ID),'ftime']
data = data[data$years<SS,]
}, error = function(e){})

LongData = data
SurvData = data.id

write.xlsx(LongData, file = paste0('longdata_s2',k,'.xlsx'))
write.xlsx(SurvData, file = paste0('survdata_s2',k,'.xlsx'))


# JM and mjoint data generation
marker.name = sort(unique(LongData$item))
WideData=list()
for (i in 1:n){
  tmp=LongData[LongData$ID==i,]
  nrm=length(unique(tmp$years))
  # long data
  wideData=matrix(0,nrm,length(marker.name))
  rest=tmp[tmp$item==marker.name[1],c(1,3:5)] # ID, years,AB.1, AB.2
  # surv data
  wideSurv=SurvData[SurvData$ID==i,]
  wideSurvData=matrix(rep(as.matrix(wideSurv),each=nrm),nrow=nrm)
  for (j in 1:length(marker.name)){
    wideData[,j]=tmp[tmp$item==marker.name[j],"value"]
  }
  colnames(wideData)=marker.name
  colnames(wideSurvData)=colnames(SurvData) # ID, fstat, ftime, x
  WideData[[i]]=cbind(rest,wideData,wideSurvData)
}
data.mjoint = data.frame(do.call(rbind,WideData))
save(data.mjoint, file = paste0('mjoint_s2',k,'.RData'))

k = k + 1
}

#-------------------------Scenario 3-------------------------#
# n = 500
# K = 3
# Sigma = AR(1)



sig = 0.1 # sd for epsilon
Ngene = 3 # number of genes
Ninfo = Ngene # signal=Ngene, no noise
eff = 1  # gene coefficients, true value of gam.y=1, gam.v=1
n = 500 # number of subjects
lam = 3 # Weibull parameter
theta = 10 # Weibull parameter
Sig = 0.5 # sd for MVN Sigma
rho = 0.5 # # paramater of AR(1)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
} # function to generate AR(1) correlation matrix

Sigma = Sig^2*ar1_cor(2, 0.5)

for (k in 1:100) {

set.seed(k) # seed setting makes the simulation data reproducible

AB1.list = lapply(1:n,function(i){
  matrix(mvrnorm(n = Ngene, rep(1.5, 2), Sigma = Sigma), nrow = Ngene)
}) # generate random intercept and slope for each gene and each subjects


Y = rep(0,n)
x = runif(n,-1,1)
beta = 1

for(i in 1:n){
  AB1 = AB1.list[[1]][1:Ninfo,]
  Atrue = sum(AB1[,1]) # summation of random intercept
  Btrue = sum(AB1[,2]) # summation of random slope
  
  logr = log(runif(1))
  ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])} # hazard function h(t)
  f = function(x){ integrate(ff,0,x)$value + logr } # Hazard function H(t)
  Y[i] = uniroot(f,c(0,50))$root
}

C = rweibull(n,2,3)
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)


len = 0.2 # the interval between each time points, also known as delta
#S = sapply(1:n,function(i){c(0,sort(runif(mmm-1,0,1))) })
data = lapply(1:n,function(i){
  
  ss=seq(0,data.id$ftime[i],by=len)
  ss = jitter(ss)
  ss[1]=0
  
  years = rep(ss,Ngene)
  item = rep(paste("gene",1:Ngene,sep=""),each=length(ss))
  ABtmp = AB1.list[[i]][rep(1:Ngene,each=length(ss)),]
  
  data.frame(ID=i,item=item,years=years,AB=ABtmp)
})

data = do.call(rbind,data)
value = data[,'AB.1'] + data[,'AB.2']*data[,'years'] + rnorm(nrow(data),0,sig)
data = cbind(data,value)

SS = data.id[match(data$ID,data.id$ID),'ftime']
data = data[data$years<SS,]


LongData = data
SurvData = data.id
write.xlsx(LongData, file = paste0('longdata_s3',k,'.xlsx'))
write.xlsx(SurvData, file = paste0('survdata_s3',k,'.xlsx'))

# JM and mjoint data generation
marker.name = sort(unique(LongData$item))
WideData=list()
for (i in 1:n){
  tmp=LongData[LongData$ID==i,]
  nrm=length(unique(tmp$years))
  # long data
  wideData=matrix(0,nrm,length(marker.name))
  rest=tmp[tmp$item==marker.name[1],c(1,3:5)] # ID, years,AB.1, AB.2
  # surv data
  wideSurv=SurvData[SurvData$ID==i,]
  wideSurvData=matrix(rep(as.matrix(wideSurv),each=nrm),nrow=nrm)
  for (j in 1:length(marker.name)){
    wideData[,j]=tmp[tmp$item==marker.name[j],"value"]
  }
  colnames(wideData)=marker.name
  colnames(wideSurvData)=colnames(SurvData) # ID, fstat, ftime, x
  WideData[[i]]=cbind(rest,wideData,wideSurvData)
}
data.mjoint = data.frame(do.call(rbind,WideData))
save(data.mjoint, file = paste0('mjoint_s3_', k, '.RData'))
}
