####Spec function####
spec_function <- function(rolist){
  if (rolist[["model_specification"]][["threshhold_included"]] == TRUE ) {
    model_oil <- "fGARCH"
    submodel_oil <- "TGARCH"
  }else{
    model_oil <- "sGARCH"
    submodel_oil <- NULL
  }
  
  order_ar <- rolist[["model_specification"]][["number_ar"]]
  order_ma <- rolist[["model_specification"]][["number_ma"]]
  
  if (rolist[["model_specification"]][["distribution"]]=="t"){
    distr_oil <- "std"
  }else{distr_oil <- "norm"}
  
  
  specs <- ugarchspec(variance.model=list(model=model_oil, submodel = submodel_oil, 
                                          garchOrder=c(order_ar,order_ma)),
                      mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                      distribution.model=distr_oil)
  specs
}



####DCC####
dcc_function <- function(rub_list, oil_list){
  rub_specs <- spec_function(rub_list)
  oil_specs <- spec_function(oil_list)
  specs <- list(rub_specs,oil_specs)
  returns <- cbind(rub_list[["return_data"]],oil_list[["return_data"]])
  uspec.n = multispec(specs)
  multf <- multifit(uspec.n, returns, solver ='solnp')
  coefs_oil <- as.vector(t(oil_list[["garch_coefs"]][1,]))
  names(coefs_oil) <- colnames(oil_list[["garch_coefs"]])
  
  coefs_rub<- as.vector(t(rub_list[["garch_coefs"]][1,]))
  names(coefs_rub) <- colnames(rub_list[["garch_coefs"]])
  
  multf@fit[[1]]@fit[["coef"]] <- coefs_rub#rub
  multf@fit[[2]]@fit[["coef"]] <- coefs_oil#oil
  
  multf@fit[[1]]@fit[["solver"]][["sol"]][["pars"]] <- coefs_rub#rub
  multf@fit[[2]]@fit[["solver"]][["sol"]][["pars"]] <- coefs_oil#oil
  
  for (i in 1:ncol(rub_list[["garch_coefs"]])) {
    multf@fit[[1]]@fit[["ipars"]][which(rownames(multf@fit[[1]]@fit[["ipars"]]) == colnames(rub_list[["garch_coefs"]])[i]),1] <- rub_list[["garch_coefs"]][1,i]
  }  
  
  for (i in 1:ncol(oil_list[["garch_coefs"]])) {
    multf@fit[[2]]@fit[["ipars"]][which(rownames(multf@fit[[2]]@fit[["ipars"]]) == colnames(oil_list[["garch_coefs"]])[i]),1] <- oil_list[["garch_coefs"]][1,i]
  }  
  
  
  for (i in 1:ncol(rub_list[["garch_coefs"]])) {
    multf@fit[[1]]@fit[["matcoef"]][which(rownames(multf@fit[[1]]@fit[["matcoef"]]) == colnames(rub_list[["garch_coefs"]])[i]),1] <- rub_list[["garch_coefs"]][1,i]
  }  
  
  for (i in 1:ncol(oil_list[["garch_coefs"]])) {
    multf@fit[[2]]@fit[["matcoef"]][which(rownames(multf@fit[[2]]@fit[["matcoef"]]) == colnames(oil_list[["garch_coefs"]])[i]),1] <- oil_list[["garch_coefs"]][1,i]
  }  
  
  
  for (i in 1:ncol(rub_list[["garch_coefs"]])) {
    multf@fit[[1]]@model[["pars"]][which(rownames(multf@fit[[1]]@model[["pars"]]) == colnames(rub_list[["garch_coefs"]])[i]),1] <- rub_list[["garch_coefs"]][1,i]
  }  
  
  for (i in 1:ncol(oil_list[["garch_coefs"]])) {
    multf@fit[[2]]@model[["pars"]][which(rownames(multf@fit[[2]]@model[["pars"]]) == colnames(oil_list[["garch_coefs"]])[i]),1] <- oil_list[["garch_coefs"]][1,i]
  }  
  
  spec1 <- dccspec(uspec = uspec.n, dccOrder = c(1,1), distribution = "mvt")#, 'mvt'fixed.pars = "fixed.se")#fixed.pars = as.list(coef(sgarch.fit)))
  fit1 <- dccfit(spec1, data = returns, solver = "nlminb", fit = multf, fit.control = list(scale =TRUE)) #solver =c("solnp", "nlminb", "lbfgs", "gosolnp"))
  fit1
}

####Call Function####
full_sample <- readRDS("C:/Users/johan/Documents/GitHub/Volatile-Gamma/output/univariateModels/univariate_garchs_full_sample.rds")
  rub_list <- full_sample[["rub"]]
  oil_list <- full_sample[["oil"]]
  
fitdcc <- dcc_function(full_sample[["rub"]],full_sample[["oil"]])
fitdcc
