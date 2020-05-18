full_sample <- readRDS("C:/Users/johan/Documents/GitHub/Volatile-Gamma/output/univariateModels/univariate_garchs_full_sample.rds")


####Oil####
if (full_sample[["oil"]][["model_specification"]][["threshhold_included"]] == TRUE ) {
  model_oil <- "fGARCH"
  submodel_oil <- "TGARCH"
}else{
  model_oil <- "sGARCH"
  submodel_oil <- NULL
  }

order_ar <- full_sample[["oil"]][["model_specification"]][["number_ar"]]
order_ma <- full_sample[["oil"]][["model_specification"]][["number_ma"]]

if (full_sample[["oil"]][["model_specification"]][["distribution"]]=="t"){
  distr_oil <- "std"
}else{distr_oil <- "norm"}


oil_spec <- ugarchspec(variance.model=list(model=model_oil, submodel = submodel_oil, 
                                           garchOrder=c(order_ar,order_ma)),
                       mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                       distribution.model=distr_oil)


returns_oil <- full_sample[["oil"]][["return_data"]]
ug_fit = ugarchfit(spec = oil_spec, data = returns_oil, solver ='hybrid')



####Rub####
if (full_sample[["rub"]][["model_specification"]][["threshhold_included"]] == TRUE ) {
  model_rub <- "fGARCH"
  submodel_rub <- "TGARCH"
}else{
  model_rub <- "sGARCH"
  submodel_rub <- NULL
}

order_ar <- full_sample[["rub"]][["model_specification"]][["number_ar"]]
order_ma <- full_sample[["rub"]][["model_specification"]][["number_ma"]]

if (full_sample[["rub"]][["model_specification"]][["distribution"]]=="t"){
  distr_rub <- "std"
}else{distr_rub <- "norm"}


rub_spec <- ugarchspec(variance.model=list(model=model_rub, submodel = submodel_rub, 
                                           garchOrder=c(order_ar,order_ma)),
                       mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                       distribution.model=distr_rub)


returns_rub <- full_sample[["rub"]][["return_data"]]
ug_fit = ugarchfit(spec = rub_spec, data = returns_rub, solver ='hybrid', fit.control = list(stationarity = 1, fixed.se = 0, scale = 0))

####DCC Prep####
specs <- list(rub_spec,oil_spec)
returns <- cbind(returns_rub,returns_oil)
uspec.n = multispec(specs)
multf <- multifit(uspec.n, returns, solver ='solnp')

fixed.pars.oil <- full_sample[["oil"]][["garch_coefs"]]
fixed.pars.rub <- full_sample[["rub"]][["garch_coefs"]]

multf@fit[[1]]@model[["pars"]][1,1] <- fixed.pars.rub$mu_return
multf@fit[[1]]@model[["pars"]][7,1] <- fixed.pars.rub$constant_garch
multf@fit[[1]]@model[["pars"]][8,1] <- fixed.pars.rub$ma1
multf@fit[[1]]@model[["pars"]][9,1] <- fixed.pars.rub$ar1
multf@fit[[1]]@model[["pars"]][17,1] <- fixed.pars.rub$df_t_distrib

multf@fit[[2]]@model[["pars"]][1,1] <- fixed.pars.oil$mu_return
multf@fit[[2]]@model[["pars"]][7,1] <- fixed.pars.oil$constant_garch
multf@fit[[2]]@model[["pars"]][8,1] <- fixed.pars.oil$ma1
multf@fit[[2]]@model[["pars"]][9,1] <- fixed.pars.oil$ar1
multf@fit[[2]]@model[["pars"]][11,1] <- fixed.pars.oil$threshhold_coef
multf@fit[[2]]@model[["pars"]][17,1] <- fixed.pars.oil$df_t_distrib

spec1 <- dccspec(uspec = uspec.n, dccOrder = c(1,1), distribution = "mvt")#, 'mvt'fixed.pars = "fixed.se")#fixed.pars = as.list(coef(sgarch.fit)))
fit1 <- dccfit(spec1, data = returns, solver = "solnp", fit = multf, fit.control = list(scale =TRUE)) #solver =c("solnp", "nlminb", "lbfgs", "gosolnp"))




