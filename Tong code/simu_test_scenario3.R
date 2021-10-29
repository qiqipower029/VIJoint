## simulate data ###
# rm(list=ls())
require(survival)
require(splines)


#sig = 0.1
Ngene = 5 # total genes; C; 3,5,10
Ninfo = Ngene # signal=Ngene, no noise
eff = 3/Ngene # gene coefficients, true value of gam.y=1, gam.v=1
n = 500 
#mmm = 5 
#nt = sample(2:mmm, n, replace=TRUE)
#ID = rep(1:n, nt)
lam = 3 # weibull distribution parameters
theta = 10

AB1.list = lapply(1:n,function(i){
  cbind(runif(Ngene,0,1), runif(Ngene,0,1))
}) # random effects of each gene

Y = rep(0,n)
x = runif(n,-1,1)
beta = 1 # gam.v

for(i in 1:n){ #survival outcomes
  AB1 = AB1.list[[i]][1:Ninfo,]
  Atrue = sum(AB1[,1]) # random intercept
  Btrue = sum(AB1[,2]) # random slope
  
  logr = log(runif(1))
  ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])}
  f = function(x){ integrate(ff,0,x)$value + logr }
  Y[i] = uniroot(f,c(0,50))$root
}

C = rweibull(n,2,3) # D, number of events can not be too small
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)# x baseline covarites

# longitudial data
len = 0.5 # time interval 
#S = sapply(1:n,function(i){c(0,sort(runif(mmm-1,0,1))) })
data = lapply(1:n,function(i){
  
  ss=seq(0,data.id$ftime[i],by=len)
  ss = jitter(ss) # longitudinal follow-up time
  ss[1]=0 
  
  years = rep(ss,Ngene)
  item = rep(paste("gene",1:Ngene,sep=""),each=length(ss))
  ABtmp = AB1.list[[i]][rep(1:Ngene,each=length(ss)),] # long-long format data
  
  data.frame(ID=i,item=item,years=years,AB=ABtmp)
})

data = do.call(rbind,data)

# use LMM to obtain longitudinal responses
# Vary variances of genes; 
sig= seq(0.1,1,.2)
item_nm=paste("gene",1:Ngene,sep="")
sel=lapply(1:length(item_nm),function(ii){
  data$item==item_nm[ii]
})
value=numeric(nrow(data))
for (i in 1:Ngene){
  nrw=length(which(sel[[i]]==TRUE))
  value[sel[[i]]] = data[sel[[i]],'AB.1'] + data[sel[[i]],'AB.2']*data[sel[[i]],'years'] + rnorm(nrw,0,sig[i])
}
#value = data[,'AB.1'] + data[,'AB.2']*data[,'years'] + rnorm(nrow(data),0,sig)
data = cbind(data,value)

SS = data.id[match(data$ID,data.id$ID),'ftime'] # force "surv data" match the order of subjects in long data
data = data[data$years<SS,]


LongData = data
SurvData = data.id