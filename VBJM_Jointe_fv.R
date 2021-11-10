#rm(list=ls())

require(Rcpp)
require(VBJM)
require(joineRML)

args=commandArgs(trailingOnly = TRUE) # Rscript i scenario1
print(args)
ffs=as.integer(args[1])
print(ffs)

source("/projects/sph_sbasu/tzhang87/JM/simu_test_scenario1.R")
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
    "gene3" = gene3 ~ years,
    "gene4" = gene4 ~ years,
    "gene5" = gene5 ~ years),
  formLongRandom = list( "gene1" = ~ years | ID, 
                         "gene2" = ~ years | ID, 
                         "gene3" = ~ years | ID,
                         "gene4" = ~ years | ID,
                         "gene5" = ~ years | ID),
  formSurv = Surv(ftime, fstat) ~ x, data = data.JM,
  timeVar = "years",
  control = list(tol0 = 0.001, burnin = 400)
))

#system.time(fit.mjoint.gn<-update(fit.mjoint,gammaOpt = "GN"))

summary(fit.mjoint)
Surv.coef.CI=confint(fit.mjoint)[11:16,]

filename=paste("/projects/sph_sbasu/tzhang87/JM/result_5gene_sce1",ffs,".rda",sep="")
save(res, res_summary, fit.mjoint, Surv.coef.CI, time_VBJM, time_JM, Delta, corv, file=filename)

quit(save="no")


#772372