# Tree Garch - Test

#### Pre
# librarys
  library(readxl)    
  library(readr)   
  library(data.table)
  library(tidyverse)  
  library(stats)
  library(ggfortify)  
  library(zoo)        
  library(xtable)     
  library(lmtest)     
  library(forecast) 
  library(tseries)
  library(xts)
  library(quantmod)
  library(fGarch)
  library(rugarch)
  library(rmgarch)
  library(psych)
  library(MASS)
  library(fitdistrplus)
  
  filter <- dplyr::filter
  select <- dplyr::select
  setwd("~/GitHub/Volatile-Gamma") # setwd
  
  source("scripts/functions.R") # functions
  source("scripts/garchFunction.R") # functions

# Data Import ----
  # load arima-errors
  ts_r = read.table('data/data_e.csv', sep = ',')
  ts_r = xts(ts_r, order.by = as.Date(rownames(ts_r)))
  
  
## remove later
  quantile_outliers = 0.002 #disabled
  ts_r <- ts_r[(ts_r$rub > quantile(ts_r$rub,quantile_outliers)) & (ts_r$rub < quantile(ts_r$rub,1-quantile_outliers)) & (ts_r$oil > quantile(ts_r$oil,quantile_outliers))& (ts_r$oil < quantile(ts_r$oil,1-quantile_outliers))]
  
          

# Garch Function ----
    # to see number of lags: just do 
          
    # ACF / PACF r^2 
    #need to adjust confidence bounds, portm-test
      ts_r2 <- ts_r^2
      acf(ts_r2$oil_errors, lag.max = 30, plot = TRUE)
      acf(ts_r2$rub_errors, lag.max = 30, plot = TRUE)  
      pacf(ts_r2$oil_errors, lag.max = 30, plot = TRUE)
      pacf(ts_r2$rub_errors, lag.max = 30, plot = TRUE)  
      # suggests high order AR lags
          
          
    # test for asymmetries
      sampleAutocorrelation(ts_r$oil_errors, "oil", 0.05)
      sampleAutocorrelation(ts_r$rub_errors, "RUBUSD", 0.05)

          
     # sign bias tests   
     # TODO
      a=(sp500<0) #Sign Bias
      a=c()
      for (i in 1:length(ibm))
      {
        a=c(a,100*min(sp500[i],0)) #Negative Size Bias
        #a=c(a,100*max(sp500[i],0)) #Positive Size Bias
      }

      h=5
      b=lm(100*sp500[(h+1):length(ibm)]^2~a[1:(length(ibm)-h)])
      summary(b)
          
          
     
   # tree garch test ----
            returns = ts_r$oil_errors

              # define possible split variables
                # past returns, epsilons, variances of own process and other processes
                 means = colMeans(ts_r)
                 epsilon = sweep(ts_r,2,means)
                 epsilon_sq = epsilon^2
                 colnames(epsilon) = paste0(colnames(ts_r),"_epsilon")
                 colnames(epsilon_sq) = paste0(colnames(ts_r),"_epsilon_sq")
                 
                 
                 base_split_variables = as.xts(cbind(ts_r, epsilon, epsilon_sq))
                 
                 # get 2 lags of each variable. In VAR tests, dependencies were not consistent above the 2nd lag. Are excluded to get parsimonious computationally feasible model
                 max_lags = 2 # number lags for loop such that for all  split vars the same dataset is used (NA in the first observations otherwise)
                 lag1 = lag(base_split_variables,1)
                 lag2 = lag(base_split_variables,2)
                 colnames(lag1) = paste0(colnames(base_split_variables),"_lag1")
                 colnames(lag2) = paste0(colnames(base_split_variables),"_lag2")
                 
                 split_variables = as.data.frame(cbind(lag1 , lag2)) # data frame instead of time series

                # possible others: interactions?
                
                 
             # step 1) find optimal GARCH 1/1 for full sample. remove first max_lags obs since they are not used by TreeGarch either
                 source("scripts/garchFunction.R") # functions
                 # inputs fct
                 returns=ts_r$rub_errors
                 ar = 1
                 ma = 1
                 threshhold = F
                 th_value  = 0 # not optimized within fct
                 data_threshhold = 0
                 type = "GARCH"
                 distribution ="normal"
                 start_parms = c(0, 0.1, 0.9,0.1)
                 opt_parms= nlm(garchEstimation,start_parms,
                                returns = returns[(max_lags+1):length(returns)],  ar = ar, ma = ma,
                                threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                                type=type, distribution=distribution,
                                print.level=2,iterlim=1000, check.analyticals=1)
               

              
          # step 2) split sample via reduction in log likelihood
               vector_quantiles = seq(1, 7)*0.125 # for threshholds
             
             # choose variable for partition
               list_split_variables = colnames(split_variables)

               split_one_covariate = drop_na(select(split_variables, list_split_variables[2])) # assign current split variable and remove first missing obs
               split_threshholds = quantile(split_one_covariate[,1],vector_quantiles)

               
             # starting values are parms of first GARCH
               start_parms = opt_parms$estimate  
              
              for (i in 1:length(split_threshholds)) {
                # split sample starting from first obs that is not NA for splitting value
                subsample1 = returns[(split_one_covariate >=split_threshholds[i])[(max_lags+1):length(returns)]] 
                subsample2 = returns[(split_one_covariate <split_threshholds[i])[(max_lags+1):length(returns)]]
                
                # estimate GARCH in subsamples and get sum of likelihood
                opt_parms1= nlm(garchEstimation,start_parms,
                               returns = subsample1,  ar = ar, ma = ma,
                               threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                               type=type, distribution=distribution,
                               iterlim=1000, check.analyticals=1)
                opt_parms2= nlm(garchEstimation,start_parms,
                                returns = subsample2,  ar = ar, ma = ma,
                                threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                                type=type, distribution=distribution,
                                iterlim=1000, check.analyticals=1)
                
                print(opt_parms$minimum)
                print(opt_parms1$minimum+ opt_parms2$minimum )
                
              }
              
            
              # step 3) prune: choose subtree that minimized sum AIC
            
           

   
          # news impact curve:
            # IMPLEMENT
               
               
          # compare results of audrino and our fct ----
          garchEstimation(theta=start_parms, returns, ar, ma, threshhold,th_value,data_threshhold,type, distribution)
          my.loglike.t( my.optpar$estimate)
          
          source("scripts/garchFunction.R") # functions
          # inputs fct
            returns=ts_r$rub_errors
            ma = 1
            ar = 2
            threshhold = T
            th_value  = 0 # not optimized within fct
            data_threshhold = 0 # not implemented 
            distribution ="t"
            
          # starting parms
            start_parms = c(0,0.1,  rep(0.1,ma), rep(0.5,ar)) # initialize parms. 
            if(threshhold==T){
              start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
            }
            if(distribution=="t"){
              start_parms=  c(start_parms, 6) # keep df_t > 2 due to likelihood fct
            }
            
          # get list for input specification
            number_parms_estimated = length(start_parms)
            model_specification = list(number_parms_estimated,ma,ar,threshhold, distribution,th_value)
            names(model_specification) = c("number_parms_estimated","number_ma","number_ar","threshhold_included", "distribution","th_value")
            

            
            opt_parms= nlm(garchEstimation,start_parms,
                           returns = returns,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                           threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                           distribution=model_specification$distribution,
                           print.level=1,steptol = 1e-6, iterlim=1000, check.analyticals=T)
          
          # model parameters
            names_parms =  c("mu_return", "constant_garch", paste0("ma",seq(1:ma)), paste0("ar",seq(1:ar)))
            if(threshhold==T){
              names_parms=  c(names_parms, "threshhold_coef") # set asymmetry parameter to 0 
            }
            if(distribution=="t"){
              names_parms=  c(names_parms, "df_t_distrib") # keep df_t > 2 due to likelihood fct
            }
            
            garch_coefs = as.data.frame(t(c(opt_parms$estimate[1:2], opt_parms$estimate[3:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma):length(opt_parms$estimate)])))
            colnames(garch_coefs) = names_parms
           
            print("GARCH coefs") 
            print(garch_coefs)
            
          # model evaluation
            
            # stationarity
              sum_coefs = sum(garch_coefs[,3:(2+ar+ma)])
              if(threshhold==T){
                sum_coefs=  sum_coefs + sum(returns<=th_value)/(length(returns))*garch_coefs$threshhold_coef # adjust if threshhold is active
              }
              print("Sum of Coefs (including threshhold)")
              print(sum_coefs)
          
            # model selection
              print("log_likelihood")
              loglik_model =-opt_parms$minimum
              loglik_model
              print("AIC")
              aic_model = my_aic(loglik_model, model_specification$number_parms_estimated)
              aic_model
              print("BIC")
              bic_model = my_bic(loglik_model, model_specification$number_parms_estimated, length(returns))
              bic_model
              
          # save model evaluation in list  
            model_evaluation = list(sum_coefs,loglik_model, aic_model, bic_model)
            names(model_evaluation) = c("sum_ar_ma_coefs","log_lik","aic_model","bic_model")
          
          # save model evaluation in list
            garch_model = list(colnames(returns), returns, garch_coefs, model_specification, model_evaluation)
            names(garch_model) = c("series_name", "return_data", "garch_coefs", "model_specification", "model_evaluation")
          
            saveRDS(garch_model, file = "output/garch_oil.rds")
            test = readRDS("output/garch_oil.rds")
            
          m1=garchFit(returns~garch(1,2),data=returns,cond.dist = "std",trace=F)
          summary(m1)
          
          
  
          
          
          ####Univariate Garch Opt function ####
      
          # ug_spec= ugarchspec(variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(0, ), include.mean = TRUE), distribution.model="norm")
          # ug_fit= ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
          # ug_fit

          opt_garch<- function(ar,ma,returns) {
            AIC <- matrix(data=NA,nrow=2,ncol=2)
            for (j in 1:1) {
              for (i in 1:1) {
                ug_spec <- ugarchspec(variance.model=list(model="fGARCH", submodel= "TGARCH",garchOrder=c(i,j)), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="std") #,submodel="TGARCH"
                ugfit = ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
                info <- infocriteria(ugfit)
                AIC[j,i] <- info[2,1]
              }
            }
            opt_garch <- which(AIC == min(AIC,na.rm=TRUE), arr.ind=TRUE)
            ug_spec <- ugarchspec(variance.model=list(model="fGARCH",submodel= "TGARCH", garchOrder=c(opt_garch[1],opt_garch[2])), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="std") # ,submodel="TGARCH"
            ug_fit = ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
            output <- list("Specs"=ug_spec,"fit"=ug_fit)
          }
          
          ####Uni-Garch Specs Oil and RUB####
          output_oil <- opt_garch(0,0, ts_r$oil_errors)
          output_rub <- opt_garch(0,0, ts_r$rub_errors)
          
          #Garch summary Oil
          output_oil[["fit"]]
          
          #Garch summary Oil
          output_rub[["fit"]]
          
          ####Bivar Garch####
          cor(ts_r$rub_errors,ts_r$oil_errors)
          var(ts_r)
          var(ts_r^2)
          cor(ts_r$rub_errors^2,ts_r$oil_errors^2)
          
