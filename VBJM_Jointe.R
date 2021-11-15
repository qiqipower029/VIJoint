#rm(list=ls())

require(Rcpp)
require(VBJM)
require(joineRML)

args=commandArgs(trailingOnly = TRUE) # Rscript i scenario1
print(args)
ffs=as.integer(args[1])
print(ffs)

source("/projects/sph_sbasu/tzhang87/JM/simu_test.R")
source("/projects/sph_sbasu/tzhang87/JM/VBJM_help.R")

#sourceCpp("VBJM.cpp")

### VBJM ###
#### initiation ####
#load("/projects/sph_sbasu/tzhang87/JM/simu1.RData")

init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5) #quadrature pts

##
time_VBJM=system.time({
  res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})

res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)


### mjoint ###

marker.name = sort(unique(LongData$item))
#data.tmp = LongData[LongData$item==marker.name[i],]
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
data.JM=data.frame(do.call(rbind,WideData))


time_JM=system.time(fit.mjoint <- mjoint(
  formLongFixed = list(
    "gene1" = gene1 ~ years,
    "gene2" = gene2 ~ years,
    "gene3" = gene3 ~ years),
  formLongRandom = list( "gene1" = ~ years | ID, 
                         "gene2" = ~ years | ID, 
                         "gene3" = ~ years | ID),
  formSurv = Surv(ftime, fstat) ~ x, data = data.JM,
  timeVar = "years",
  control = list(tol0 = 0.001, burnin = 400)
))

#system.time(fit.mjoint.gn<-update(fit.mjoint,gammaOpt = "GN"))

summary(fit.mjoint)
Surv.coef.CI=confint(fit.mjoint)[7:10,]

# time=rbind(time_VBJM,time_JM)
# rownames(time)=c("VBJM","MCJM")
# gammas=fit.mjoint$coefficients$gamma
# gam.sd=summary(fit.mjoint)[[2]][,2]
# fit.mjoint$Hessian


filename=paste("/projects/sph_sbasu/tzhang87/JM/result_3gene",ffs,".rda",sep="")
save(res, res_summary, fit.mjoint, Surv.coef.CI, time_VBJM, time_JM, Delta, file=filename)

quit(save="no")


## Wrap up results
# result_sum=function(ngene,sce){
#   iter=100
#   tm_JM=tm_VBJM=matrix(0,iter,5)
#   tab_Delta=matrix(0,iter,2)
#   gamma1=numeric(iter)
#   JM_CI_res=matrix(0,iter, 1+ngene) # 1+5, 1+10
#   VBJM_CI_res=matrix(0,iter, 1+ngene) # 1+5, 1+10
#   VBJM_CI=JM_CI=vector(mode="list",length = 100)  # <------
#   #VBJM_CI=JM_CI=list()
#   # true x always = 1; 3gene true gamma=1; 5gene true gamma=3/5; 10gene true gamma 3/10;
#   true_x=1
#   true_gamma=c(1,3/5,.3)
#   true_gam=true_gamma[sce] # sce=gene3,gene5, gene10
#   locParm.VBJM=list(c(7:10), c(11:16), c(21:31))
#   error=NULL
# 
#   for (i in 1:100){
#     #tryCatch({
#     #filename=paste("VBJM_res/res_",ngene,"gene/result_",ngene,"gene", i, ".rda",sep="") # 5gene, 10gene
#     filename=paste("VBJM_res/res_",ngene,"gene_sce3/result_",ngene,"gene_sce3", i, ".rda",sep="") # 5gene, 10gene
#     load(filename)
#     tm_JM[i,]=time_JM
#     colnames(tm_JM)=names(time_JM)
#     tm_VBJM[i,]=time_VBJM
#     colnames(tm_VBJM)=names(time_VBJM)
#     tab_Delta[i,]=table(Delta)
#     colnames(tab_Delta)=names(table(Delta))
#     # check JM CI coverage
#     colnames(JM_CI_res)=rownames(Surv.coef.CI)
#     for (j in 1:nrow(Surv.coef.CI)){
#       if (j==1){true=true_x} else{true=true_gam}
#       if (Surv.coef.CI[j,1]<= true & true<=Surv.coef.CI[j,2]){JM_CI_res[i,j]=1} else {JM_CI_res[i,j]=0}
#     }
#     # check VBJM CI coverage
#     locParm=locParm.VBJM[[sce]]
#     colnames(VBJM_CI_res)=rownames(res_summary)[locParm]
#     for (y in locParm){
#       z=which(locParm==y)
#       if (z==1){true=true_x} else{true=true_gam}
#       if (res_summary[y,3]<= true & true<=res_summary[y,4]){VBJM_CI_res[i,z]=1} else {VBJM_CI_res[i,z]=0}
#     }
#     VBJM_CI[[i]]=res_summary[locParm,3:4]
#     JM_CI[[i]]=Surv.coef.CI
#     #}, warning=function(w){ print(w)})
#     
#     sd=solve(fit.mjoint$Hessian)
#     sigmak=sqrt(diag(sd[66:70,66:70]))
#   }
#   CI_vbjm=do.call(cbind,VBJM_CI)
#   CI_jm=do.call(cbind,JM_CI)
#   # return(list(head(JM_CI_res),head(CI_jm),head(VBJM_CI_res),head(CI_vbjm),
#   #             head(tm_JM),head(tm_VBJM), head(tab_Delta)))
# 
#   return(list(JM_CI_res=JM_CI_res,VBJM_CI_res=VBJM_CI_res, 
#           JM_Time=tm_JM,VBJM_Time=tm_VBJM, DeltaFreq=tab_Delta))
# }
# 
# res3=result_sum(3,1)
# JM3CI=apply(res3$JM_CI_res,2,sum)
# VB3JMCI=apply(res3$VBJM_CI_res,2,sum)
# tm3JM=apply(res3$JM_Time,2,mean) # check zero
# tm3VBJM=apply(res3$VBJM_Time,2,mean)
# apply(res3$DeltaFreq,2,median)
# apply(res3$DeltaFreq,2,min)
# 
# res5=result_sum(5,2)
# JM5CI=apply(res5$JM_CI_res,2,sum)
# VB5JMCI=apply(res5$VBJM_CI_res,2,sum)
# tm5JM=apply(res5$JM_Time,2,mean) # check zero
# tm5VBJM=apply(res5$VBJM_Time,2,mean)
# apply(res5$DeltaFreq,2,median)
# apply(res5$DeltaFreq,2,min)
# 
# MIS=c(26,28:36, 38:40, 45, 48, 51, 53, 55, 56, 59, 65, 70, 74:76, 78, 81, 84, 89:91, 94:96, 99)
# vec=1:135
# nMIS=vec[-MIS]
# 
# res10=result_sum(10,3)
# JM10CI=apply(JM_CI_res,2,sum)
# VB10JMCI=apply(VBJM_CI_res,2,sum)
# tm10JM=apply(tm_JM,2,mean) # check zero
# tm10VBJM=apply(tm_VBJM,2,mean)
# apply(tab_Delta,2,median)
# apply(tab_Delta,2,min)

# JM3CI=c(90,92,94,93)
# VBJM3CI=c(96,97,98,99)
# 
# plot(1:11, 90:100, main="Numbers of CI covers the true value")
# lines(1:4, JM3CI, col="red", lty=2)
# lines(1:4, VBJM3CI, col="red", lty=1)
# nm.jm=rownames(Surv.coef.CI)
# sd=solve(fit.mjoint$Hessian)
# sigmak_sd=sqrt(diag(sd[66:70,66:70]))
# sigmak=sqrt(fit.mjoint$coefficients$sigma2)
# SurvcoefCI=cbind(sigmak+sigmak_sd*qnorm(.025),sigmak+sigmak_sd*qnorm(.975))
# Surv.coef.CI=rbind(Surv.coef.CI,SurvcoefCI)
# rownames(Surv.coef.CI)=c(nm.jm,paste("sigmak",1:ngene,sep=""))
# # nm.vbjm=rownames(res_summary)[locParm]
# # SD=-1*solve(res$H)
# # sig2k=sqrt(res$sig2)
# # sig2k_sd=
# # res_summary[(ngene*2+3+ngene+1):(ngene*2+3+ngene*2),3]=sig2k+sig2k_sd*qnorm(.025)
# # res_summary[(ngene*2+3+ngene+1):(ngene*2+3+ngene*2),4]=sig2k+sig2k_sd*qnorm(.975)
# # rownames(res_summary)=c(nm.vbjm,paste("sig2",1:ngene,sep=""))