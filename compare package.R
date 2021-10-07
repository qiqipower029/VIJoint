#------------------------------Compare packages------------------------------#
#--------------------------------Jieqi Tu------------------------------------#

library(JM)
library(joineRML)
library(VBJM)
require(Matrix)
require(statmod)
library(tidyverse)
require(survival)
data("pbc2")
#-------------------Using simu1.RData-----------------#
#--------VBJM package---------#
load("./simu1.RData")
simudata = full_join(LongData, SurvData, by="ID") %>% as.data.frame()

init_list = VBJM_init(LongData=LongData, SurvData=SurvData, n_points=5)

##
system.time({
  res = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})


res_summary = VBJM_get_summary(marker.name =init_list$marker.name, res=res)

#--------joineRML package---------#

# Subset the data into two categories: gene1 and gene2
gene1 = simudata %>% filter(item == "gene1")
gene2 = simudata %>% filter(item == "gene2")

gene1.fit = mjoint(
  formLongFixed = list("year" = value ~ years),
  formLongRandom = list("year" = ~ years|ID),
  formSurv = Surv(ftime, fstat) ~ x,
  data = gene1,
  timeVar = "years"
)
summary(gene1.fit)

gene2.fit = mjoint(
  formLongFixed = list("year" = value ~ years),
  formLongRandom = list("year" = ~ years|ID),
  formSurv = Surv(ftime, fstat) ~ x,
  data = gene2,
  timeVar = "years"
)
summary(gene2.fit)
gene1.fit$coefficients$gamma
gene2.fit$coefficients$gamma

# 1) Gamma is very different. 2) How to check alpha using joineRML package? They may still be gamma. 
# 3) The meaning of alpha in VBJM might be misleading.


# Use wide format run mjoint()
#-------------------Using pbc2----------------#
#--------joineRML package---------#
placebo <- subset(pbc2, drug == "placebo") 
fit.pbc <- mjoint(
  formLongFixed = list(
    "bil" = log(serBilir) ~ year, "alb" = albumin ~ year, "pro" = (0.1 * prothrombin)^-4 ~ year),
  formLongRandom = list(
    "bil" = ~ year | id, "alb" = ~ year | id, "pro" = ~ year | id),
  formSurv = Surv(years, status2) ~ age,
  data = placebo,
  timeVar = "year",
  control = list(tol0 = 0.001, burnin = 400)
)
summary(fit.pbc)

#--------VBJM package---------#
long = placebo %>% select(id, year, serBilir, albumin, prothrombin)


# Transform simudata to run joineRML package
library(fastDummies)
simudata2 = simudata %>% 
  dummy_cols(select_columns = "item") %>% as.data.frame()

simu.fit = mjoint(
  formLongFixed = list("gene1" = value ~ years + item_gene1,
                       "gene2" = value ~ years
                       ),
  formLongRandom = list("gene1" = ~ item_gene1 | ID,
                        "gene2" = ~ item_gene2 | ID),
  formSurv = Surv(ftime, fstat) ~ x,
  data = simudata2,
  timeVar = "years"
)

resummary(simu.fit)
