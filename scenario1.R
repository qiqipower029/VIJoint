#------------------------Testing for Single Gene------------------------#
#------------------------Jieqi Tu------------------------#
require(VBJM)
require(JM)
require(survival)
require(Matrix)
require(statmod)
require(joineRML)
require(openxlsx)

# ----------- Create a simulation function--------------# 


VBJM_sim = function(n.sim, n) {
  est_gamma = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim)
  )
  
  est_alpha = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim)
  )
  for (i in 1:n.sim) {
    LongData = readxl::read_excel(path = paste0('longdata_s1', i, '.xlsx'))
    SurvData = readxl::read_excel(path = paste0('survdata_s1', i, '.xlsx'))
    init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5)
    res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
    res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)
    est_m[i,1] = res_summary$Estimate[[3]] 
  }
}