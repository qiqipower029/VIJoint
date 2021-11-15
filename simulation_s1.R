#------------------------Testing for Single Gene------------------------#
#------------------------Jieqi Tu------------------------#
require(VBJM)
require(JM)
require(survival)
require(Matrix)
require(statmod)
require(joineRML)

#----------VBJM package----------#
load("Scenario testing data/simu_s1_500.RData")
load("./Scenario testing data/mjoint_s1_500.RData")

init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5)


time_VBJM = system.time({
  res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})


res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)
res_summary


#----------joineRML package----------#
time_mjoint=system.time(fit.mjoint <- mjoint(
  formLongFixed = list(
    "gene1" = gene1 ~ years),
  formLongRandom = list( "gene1" = ~ years | ID),
  formSurv = Surv(ftime, fstat) ~ x, data = data.JM,
  timeVar = "years",
  control = list(tol0 = 0.001, burnin = 400)
))

summary(fit.mjoint)


#----------JM package----------#
fitLME = lme(value ~ years ,  random = ~ years | ID,
             data = LongData, control=lmeControl(opt='optim'))
fitSURV = survreg(Surv(ftime, fstat) ~  x, data = SurvData, x = TRUE)
time_JM = system.time({fitJOINT = jointModel(fitLME, fitSURV, timeVar = "years")})
summary(fitJOINT)
theta = exp(fitSURV$coefficients[1]); theta
lambda = 1/fitSURV$scale; lambda
gamma = -fitSURV$coefficients[2]/fitSURV$scale; gamma
