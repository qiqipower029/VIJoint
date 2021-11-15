## test VBJM method ##
rm(list=ls())

library(VBJM)
source("VBJM_help.R")
#require(Rcpp)
#sourceCpp("VBJM.cpp")
#### initiation ####
load("simu1.RData")

init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5)

##
system.time({
    res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})


res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)
