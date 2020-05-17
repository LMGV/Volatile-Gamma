#AR1 oil
#AR1 rub

#How does the test work
#Likelihood std for garch(2,2)
#Get the prediction funtion for dcc
  #conditional var for the forecasts
#Var

#Output for the treegarch:
  #Splits
  #list with the estimated parameters
#fixed.pars = list()) list with fixed parameters
#Rub
oil_mu      <- 0.000062    
oil_omega   <- 0.000000  
oil_alpha1  <- 0.115567    
oil_beta1   <- 0.398437    
oil_beta2   <- 0.484810    
oil_skew    <- 1.054032    
oil_shape   <- 4.534280  
fixed.pars.rub <- list(mu= oil_mu, omega = oil_omega, alpha1 = oil_alpha1, beta1 = oil_beta1,) #flexibel list
#Oil
mu      <- 0.000266    
omega   <- 0.000002    
alpha1  <- 0.044889    
beta1   <- 0.954111    
skew    <- 0.919689    
shape   <- 5.813160

#Joint
dcca1  <- 0.033145    
dccb1  <- 0.000000    

#Garch summary Oil
specs_oil <- output_oil[["Specs"]]
specs_oil@model$fixed.pars <- list(mu= oil_mu, omega = oil_omega, alpha1 = oil_alpha1, beta1 = oil_beta1...)#here add all variables

#Garch summary Oil
specs_rub <- output_rub[["Specs"]]

specs <- list(specs_rub,specs_oil)

uspec.n = multispec(specs)
multf <- multifit(uspec.n, ts_r)
multf

ctrl = list(tol = 1e-17, delta = 1e-10)
spec1 <- dccspec(uspec = uspec.n, dccOrder = c(1,1), distribution = "mvnorm")#, 'mvt'fixed.pars = "fixed.se")#fixed.pars = as.list(coef(sgarch.fit)))
fit1 <- dccfit(spec1, data = ts_r, solver = "solnp", fit = multf, solver.control = ctrl, fit.control = list(scale =TRUE)) #solver =c("solnp", "nlminb", "lbfgs", "gosolnp"))

fit1
