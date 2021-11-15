## simulate data ###
rm(list=ls())
require(survival)
require(splines)


sig = 0.1
Ngene = 2
Ninfo = 2
eff = 1
n = 200
#mmm = 5 
#nt = sample(2:mmm, n, replace=TRUE)
#ID = rep(1:n, nt)
lam = 3
theta = 10

AB1.list = lapply(1:n,function(i){
    cbind(runif(Ngene,0,1), runif(Ngene,0,1))
})

Y = rep(0,n)
x = runif(n,-1,1)
beta = 1

for(i in 1:n){
    AB1 = AB1.list[[i]][1:Ninfo,]
    Atrue = sum(AB1[,1])
    Btrue = sum(AB1[,2])
    
    logr = log(runif(1))
    ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])}
    f = function(x){ integrate(ff,0,x)$value + logr }
    Y[i] = uniroot(f,c(0,50))$root
}

C = rweibull(n,2,3)
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)


len = 0.5
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
save(LongData, SurvData, 
     file='simu1.RData')


# LongData.train = data[data$ID <= n/2,]
# LongData.test = data[data$ID > n/2,]
# 
# SurvData.train = data.id[data.id$ID<=n/2, ]
# SurvData.test = data.id[data.id$ID>n/2, ]

# save(LongData.train,LongData.test ,  
#      SurvData.train, SurvData.test, 
#      file='simu1.RData')
