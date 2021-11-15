require(survival)
require(JM)
require(Matrix)
require(statmod)

VBJM_init <- function(LongData=NULL, SurvData=NULL,n_points = 5){
    
    data.list = list()
    para.list = list()
    
    ## run JM for each biomarker to initialize 
    marker.name = sort(unique(LongData$item))
    beta.list = list()
    mu.list = list()
    V.list = list()
    sig.vec = rep(NA, length(marker.name))
    alpha.vec = rep(NA,length(marker.name))
    
    for(i in seq_along(marker.name)){
        #i = 1
        print(i)
        data.tmp = LongData[LongData$item==marker.name[i],]
        
        fitLME = lme(value ~ years ,  random = ~ years | ID,
                     data = data.tmp, control=lmeControl(opt='optim'))
        beta.list[[i]] = matrix(fitLME$coefficients$fixed,ncol=1)
        mu.list[[i]] = fitLME$coefficients$random$ID
        sig.vec[i] = fitLME$sigma
        V.list[[i]] = as.matrix(getVarCov(fitLME))
        ## JM package require to use lme ##
        ## otherwise, can use lmer as optional ##
        # fitLME = lmer(value ~ years +  (years | ID), 
        #              data = data.tmp)
        # beta.list[[i]] = matrix(summary(fitLME)$coef[,1],ncol=1)
        # sig.vec[i] = summary(fitLME)$sigma
        # mu.list[[i]] = ranef(fitLME)$ID
        # V.list[[i]] = as.matrix(bdiag(VarCorr(fitLME)))
        
        fitSURV = survreg(Surv(ftime, fstat) ~  x, data = SurvData, x = TRUE)
        # There are multiple ways to parameterize a Weibull distribution. The survreg 
        # function embeds it in a general location-scale family, which is a 
        # different parameterization than the rweibull function, and often leads
        # to confusion.
        #   survreg's scale  =    1/(rweibull shape)
        #   survreg's intercept = log(rweibull scale)
        theta = exp(fitSURV$coefficients[1])
        lambda = 1/fitSURV$scale
        gamma = -fitSURV$coefficients[2]/fitSURV$scale
        ## for gamma; can also use coxph
        # fitSurv = coxph(Surv(ftime, fstat)~x, data=SurvData)
        # gamma = matrix(fitSurv$coefficients,ncol=1)
        
        ## joint model 
        fitJOINT = jointModel(fitLME, fitSURV, timeVar = "years") 
        alpha.vec[i] = fitJOINT$coefficients$alpha
    }
    
    V.list = lapply(V.list, function(x){
        x[1:nrow(x),1:ncol(x),drop=FALSE]
    })
    
    Sigma = as.matrix(bdiag(V.list))
    V.list = matrix(rep(list(Sigma),nrow(SurvData)))
    
    mu.list = lapply(mu.list, function(x){
        lapply(1:nrow(x),function(i){
            matrix(x[1,],ncol=1)
        })
    })
    mu.list = do.call(cbind, mu.list)
    
    para.list[["mu"]] = mu.list
    para.list[["V"]] = V.list
    para.list[["Sigma"]] = Sigma
    para.list[["sig2"]] = sig.vec
    
    para.list[["beta"]] = beta.list
    para.list[["weib"]] = c(lambda, theta)
    para.list[["gamma"]] = gamma
    para.list[["alpha"]] = alpha.vec
    
    
    ## prepare the data for running the algorithm
    
    ## Y_i
    Y.list = lapply(unique(LongData$ID), function(i){
        data.tmp = LongData[LongData$ID==i,]
        lapply(marker.name,function(x){
            matrix(data.tmp$value[data.tmp$item==x],ncol=1)
        })
    })
    
    Y.list = do.call(rbind, Y.list)
    
    ## X_i
    X.list = lapply(unique(LongData$ID), function(i){
        data.tmp = LongData[LongData$ID==i,]
        lapply(marker.name,function(x){
            cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    
    X.list = do.call(rbind, X.list)
    Z.list = X.list
    
    
    ## X_i(T_i)
    X_T.list = lapply(unique(LongData$ID), function(i){
        T_i = SurvData$ftime[SurvData$ID==i]
        lapply(marker.name,function(x){
            matrix(c(1,T_i),ncol=1)
        })
    })
    X_T.list = do.call(rbind, X_T.list)
    Z_T.list = X_T.list
    
    W = matrix(SurvData$x,ncol=1)
    
    ## X_i(t) 
    ## this depends on the number of legendre Gaussian quadrature points
    Gauss.point  = gauss.quad(n_points)
    # \int_0^{T_i} f(t)dt
    # t_node = Gauss.point$nodes *(Ti/2) + Ti/2
    # w_node = Gauss.point$weights
    # Ti/2 * sum(w_node * f(t_node))
    
    X_t.list = lapply(unique(LongData$ID), function(i){
        Ti = SurvData$ftime[SurvData$ID==i]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        lapply(marker.name,function(x){
            cbind(1, t_node)
        })
    })
    
    X_t.list = do.call(rbind, X_t.list)
    Z_t.list = X_t.list
    
    w_node.list = lapply(unique(LongData$ID), function(i){
        Ti = SurvData$ftime[SurvData$ID==i]
        Gauss.point$weights*Ti/2
    })
    t_node.list = lapply(unique(LongData$ID), function(i){
        Ti = SurvData$ftime[SurvData$ID==i]
        Gauss.point$nodes *(Ti/2) + Ti/2
    })
    
    data.list[["Y"]] = Y.list
    data.list[["X"]] = X.list
    data.list[["X_t"]] = X_t.list
    data.list[["X_T"]] = X_T.list
    data.list[["Z"]] = Z.list
    data.list[["Z_t"]] = Z_t.list
    data.list[["Z_T"]] = Z_T.list
    data.list[["W"]] = W
    data.list[["GQ_w"]] = w_node.list
    data.list[["GQ_t"]] = t_node.list
    data.list[["ftime"]] = SurvData$ftime
    data.list[["fstat"]] = SurvData$fstat
    
    list(data.list=data.list, para.list=para.list,
         marker.name=marker.name)
}


VBJM_get_summary <- function(marker.name=NULL, res=NULL){
    
    
    beta.list = lapply(seq_along(marker.name), function(i){
        beta = as.numeric(res$beta[[i]])
        coef_name = paste(marker.name[i],"_beta",1:length(beta),sep="")
        names(beta) = coef_name
        beta
    })
    
    gamma = as.numeric(res$gamma)
    names(gamma) = paste("Surv_gamma",1:length(gamma),sep="")
    
    alpha = as.numeric(res$alpha)
    names(alpha) = paste(marker.name,"_alpha", sep="")
    
    weib = as.numeric(res$weib)
    names(weib) = c("Weibull_shape","Weibull_scale")
    
    para =c(do.call(c, beta.list),  gamma, alpha, weib)
    
    cov = -1*solve(res$H)
    se = round(sqrt(diag(cov)),4)[1:length(para)]
    
    res_summary = data.frame(Estimate=para, SE=se,
                             para-1.96*se, para+1.96*se)
    colnames(res_summary)[3:4] = c("95%CI_lower","95%CI_upper")
    res_summary
}

