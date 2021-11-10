#------------------------simulation function------------------------#
#------------------------Jieqi Tu------------------------#
require(VBJM)
require(JM)
require(survival)
require(Matrix)
require(statmod)
require(joineRML)
require(openxlsx)

source(file = "VBJM_help.R")

# ----------- Create simulation functions--------------# 
VBJM_sim = function(n.sim, scenario) {
  est_gamma = data.frame(
    estimate = numeric(n.sim),
    se = numeric(n.sim),
    ci_low = numeric(n.sim),
    ci_high = numeric(n.sim),
    cover = numeric(n.sim)
  )
  
  est_alpha1 = data.frame(
    estimate = numeric(n.sim),
    se = numeric(n.sim),
    ci_low = numeric(n.sim),
    ci_high = numeric(n.sim),
    cover = numeric(n.sim)
  )
  
  
  for (i in 1:n.sim) {
    LongData = readxl::read_excel(path = paste0('longdata_s', scenario, i, '.xlsx'))
    SurvData = readxl::read_excel(path = paste0('survdata_s', scenario, i, '.xlsx'))
    init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5)
    time = system.time({
    res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
    })
    res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)
    
    # store the results for each iteration
    for (j in 1:4) {
      est_gamma[i, j] = res_summary['Surv_gamma1', j]
      est_alpha1[i, j] = res_summary['gene1_alpha', j]
    }
    
    
    # calculate coverage rate of the confidence interval
    if ((est_gamma[i, 3] <= 1) && (est_gamma[i, 4] >= 1)) {
      est_gamma[i, 5] = 1
    } else est_gamma[i, 5] = 0
    
    if ((est_alpha1[i, 3] <= 1) && (est_alpha1[i, 4] >= 1)) {
      est_alpha1[i, 5] = 1
    } else est_alpha1[i, 5] = 0
  }
  
  # write the result
  VBJM_result = data.frame(
    gamma_mean = mean(est_gamma[,1]),
    alpha1_mean = mean(est_alpha1[,1]),
    gamma_sd = sd(est_gamma[,1]),
    alpha1_sd = sd(est_alpha1[,1]),
    gamma_MESE = mean(est_gamma[,2]),
    alpha1_MESE = mean(est_alpha1[,2]),
    gamma_RMSE = sqrt(mean((est_gamma[,1]-1)^2)),
    alpha1_RMSE = sqrt(mean((est_alpha1[,1]-1)^2)),
    gamma_cp = mean(est_gamma[,5])*100,
    alpha1_cp = mean(est_alpha1[,5])*100,
    time_used = mean(time[2])
  )
  
  write.xlsx(VBJM_result, file = "VBJM.xlsx")
}

JM_sim = function(n.sim, scenario) {

  # create data frames to store results
  est_gamma = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim),
    cover = rep(NA, n.sim)
  )
  
  est_alpha1 = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim),
    cover = rep(NA, n.sim)
  )
  
  for (i in 1:n.sim) {
    LongData = readxl::read_excel(path = paste0('longdata_s', scenario, i, '.xlsx'))
    SurvData = readxl::read_excel(path = paste0('survdata_s', scenario, i, '.xlsx'))
    
    fitLME = lme(value ~ years ,  random = ~ years | ID,
                 data = LongData, control=lmeControl(opt='optim'))
    fitSURV = survreg(Surv(ftime, fstat) ~  x, data = SurvData, x = TRUE)
    
    time = system.time({
    fitJOINT = jointModel(fitLME, fitSURV, timeVar = "years") 
    })
    res = summary(fitJOINT)
    
    for (j in 1:2) {
      est_alpha1[i, j] = res$`CoefTable-Event`['Assoct', j]
      est_gamma[i, j] = res$`CoefTable-Event`['x', j]
    }
    
    est_alpha1[i,3] = est_alpha1[i,1] - qnorm(0.975)*est_alpha1[i,2]
    est_alpha1[i,4] = est_alpha1[i,1] + qnorm(0.975)*est_alpha1[i,2]
    est_gamma[i,3] = est_gamma[i,1] - qnorm(0.975)*est_gamma[i,2]
    est_gamma[i,4] = est_gamma[i,1] + qnorm(0.975)*est_gamma[i,2]
    
    # calculate coverage rate of the confidence interval
    if ((est_gamma[i, 3] <= 1) && (est_gamma[i, 4] >= 1)) {
      est_gamma[i, 5] = 1
    } else est_gamma[i, 5] = 0
    
    if ((est_alpha1[i, 3] <= 1) && (est_alpha1[i, 4] >= 1)) {
      est_alpha1[i, 5] = 1
    } else est_alpha1[i, 5] = 0
    
  }
  
  # write the result
  JM_result = data.frame(
    gamma_mean = mean(est_gamma[,1]),
    alpha1_mean = mean(est_alpha1[,1]),
    gamma_sd = sd(est_gamma[,1]),
    alpha1_sd = sd(est_alpha1[,1]),
    gamma_MESE = mean(est_gamma[,2]),
    alpha1_MESE = mean(est_alpha1[,2]),
    gamma_RMSE = sqrt(mean((est_gamma[,1]-1)^2)),
    alpha1_RMSE = sqrt(mean((est_alpha1[,1]-1)^2)),
    gamma_cp = mean(est_gamma[,5])*100,
    alpha1_cp = mean(est_alpha1[,5])*100,
    time_used = mean(time[2])
  )
  
  write.xlsx(JM_result, file = "JM.xlsx")
}


joineRML_sim = function(n.sim, scenario) {
  
  est_gamma = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim),
    cover = rep(NA, n.sim)
  )
  
  est_alpha1 = data.frame(
    estimate = rep(NA, n.sim),
    se = rep(NA, n.sim),
    ci_low = rep(NA, n.sim),
    ci_high = rep(NA, n.sim),
    cover = rep(NA, n.sim)
  )
  
  
  for (i in 1:n.sim) {
    # read in data
    load(paste0('mjoint_s', scenario, i, '.RData'))
    
    time = system.time({
    simu.fit = mjoint(
      formLongFixed = list(
        "gene1" = gene1 ~ years
      ),
      formLongRandom = list(
        "gene1" = ~ years|ID
      ),
      formSurv = Surv(ftime, fstat) ~ x,
      data = data.mjoint,
      timeVar = "years"
    )})
    
    res = summary(simu.fit)
    
    for (j in 1:2) {
      est_alpha1[i, j] = res$coefs.surv['gamma_1', j]
      est_gamma[i, j] = res$coefs.surv['x', j]
    }
    
    est_alpha1[i,3] = est_alpha1[i,1] - qnorm(0.975)*est_alpha1[i,2]
    est_alpha1[i,4] = est_alpha1[i,1] + qnorm(0.975)*est_alpha1[i,2]
    est_gamma[i,3] = est_gamma[i,1] - qnorm(0.975)*est_gamma[i,2]
    est_gamma[i,4] = est_gamma[i,1] + qnorm(0.975)*est_gamma[i,2]
    
    # calculate coverage rate of the confidence interval
    if ((est_gamma[i, 3] <= 1) && (est_gamma[i, 4] >= 1)) {
      est_gamma[i, 5] = 1
    } else est_gamma[i, 5] = 0
    
    if ((est_alpha1[i, 3] <= 1) && (est_alpha1[i, 4] >= 1)) {
      est_alpha1[i, 5] = 1
    } else est_alpha1[i, 5] = 0
  
  }
  
  # write the result
  joineRML_result = data.frame(
    gamma_mean = mean(est_gamma[,1]),
    alpha1_mean = mean(est_alpha1[,1]),
    gamma_sd = sd(est_gamma[,1]),
    alpha1_sd = sd(est_alpha1[,1]),
    gamma_MESE = mean(est_gamma[,2]),
    alpha1_MESE = mean(est_alpha1[,2]),
    gamma_RMSE = sqrt(mean((est_gamma[,1]-1)^2)),
    alpha1_RMSE = sqrt(mean((est_alpha1[,1]-1)^2)),
    gamma_cp = mean(est_gamma[,5])*100,
    alpha1_cp = mean(est_alpha1[,5])*100,
    time_used = mean(time[2])
  )
  
  write.xlsx(JM_result, file = "joineRML.xlsx")
  
}


VBJM_sim(n.sim = 3, scenario = 1)
JM_sim(n.sim = 100, scenario = 1)
joineRML_sim(n.sim = 100, scenario = 1)