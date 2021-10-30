rm(list=ls())
require(survival)
require(splines)
require(survival)
require(JM)
require(Matrix)
require(statmod)
require(MASS)

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

Sigma = matrix(c(
  0.01,0,
  0,0.01
), 2, 2) # variance-covariance matrix of betas for each subject

set.seed(2021) # seed setting makes the simulation data reproducible

AB1.list = lapply(1:n,function(i){
  matrix(mvrnorm(n = Ngene, rep(2, 2), Sigma = Sigma), nrow = 1)
})

AB1.list = lapply(1:n,function(i){
  cbind(runif(Ngene,0,1), runif(Ngene,0,1))
}) 

Y = rep(0,n)
x = runif(n,-1,1)
beta = 1

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


LongData_s1 = data
SurvData_s1 = data.id
save(LongData, SurvData, 
     file='simu_s1.RData')

# JM data generation
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
data.JM_s1 = data.frame(do.call(rbind,WideData))
