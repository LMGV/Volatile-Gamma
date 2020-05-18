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
  outpathDescriptive = "output/univariateDescriptives/"
  outpathModels =  "output/univariateModels/"
  
# Data Import ----
  # load arima-errors
  # ts_r = read.table('data/data_e.csv', sep = ',')
  ts_r = read.table('data/data_removed_outliers_1.csv', sep = ',')
  ts_r = xts(ts_r, order.by = as.Date(rownames(ts_r)))

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
      sampleAutocorrelation(ts_r$oil_errors, "oil", 0.05, outpath= outpathDescriptive)
      sampleAutocorrelation(ts_r$rub_errors, "RUBUSD", 0.05, outpath= outpathDescriptive)

          
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
               
      # choose optimal GARCH model
         # possible model specifications
         ma_choices = seq(1:3)
         ar_choices = seq(1:3)
         threshhold_choices = c(T,F)
         th_value  = 0 # not optimized within fct
         data_threshhold = 0 # not implemented 
         distribution_choices =c("normal","t")

        # input data
          returns_list=list(ts_r$rub_errors, ts_r$oil_errors)
          names(returns_list) = c("rub","oil")

        # initialize selected model list for all univeriate series
          all_selected_model = vector("list", length = length(returns_list))
          names(all_selected_model) = names(returns_list)
          
        # estimate model for all univariate inputs
          for(data_iter in 1:length(returns_list)) {
            returns = returns_list[[data_iter]]
            
            # inialize result list for one univariate series
            garch_model_selection = list()
          
          # estimate all models given choices above
            for (dist_iter in 1:length(distribution_choices)){
              distribution =distribution_choices[dist_iter] # set distribution
              
              for (th_iter in 1:length(threshhold_choices)){
                threshhold = threshhold_choices[th_iter]
                
                for(ar_iter in 1:length(ar_choices)){
                  ar = ar_choices[ar_iter]
                  
                  for(ma_iter in 1:length(ma_choices)){
                    ma = ma_choices[ma_iter]
                    
                    # GARCH model Estimation
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
                          
                        # estimate GARCH model for given specification (minimize negative loglikelihood)
                          opt_parms= nlm(garchEstimation,start_parms,
                                         returns = returns,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                                         threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                                         distribution=model_specification$distribution,
                                         print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
                        
                        # get model parameters
                          names_parms =  c("mu_return", "constant_garch", paste0("ma",seq(1:ma)), paste0("ar",seq(1:ar)))
                          if(threshhold==T){
                            names_parms=  c(names_parms, "threshhold_coef") # set asymmetry parameter to 0 
                          }
                          if(distribution=="t"){
                            names_parms=  c(names_parms, "df_t_distrib") # keep df_t > 2 due to likelihood fct
                          }
                          
                          garch_coefs = as.data.frame(t(c(opt_parms$estimate[1], opt_parms$estimate[2:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma):length(opt_parms$estimate)])))
                          colnames(garch_coefs) = names_parms
                        
                        # model evaluation
                            # stationarity
                            sum_coefs = sum(garch_coefs[,3:(2+ar+ma)])
                            if(threshhold==T){
                              sum_coefs=  sum_coefs + sum(returns<=th_value)/(length(returns))*garch_coefs$threshhold_coef # adjust if threshhold is active
                            }
                            
                            # model selection
                            loglik_model =-opt_parms$minimum
                            aic_model = my_aic(loglik_model, model_specification$number_parms_estimated)
                            bic_model = my_bic(loglik_model, model_specification$number_parms_estimated, length(returns))
                            
                            # save model evaluation in list  
                            model_evaluation = list(sum_coefs,loglik_model, aic_model, bic_model)
                            names(model_evaluation) = c("sum_ar_ma_coefs","log_lik","aic_model","bic_model")
                        
                       # add model to model selection list
                        garch_model = list(colnames(returns), returns, garch_coefs, model_specification, model_evaluation)
                        names(garch_model) = c("series_name", "return_data", "garch_coefs", "model_specification", "model_evaluation")
                        garch_model_selection[[length(garch_model_selection)+1]] = garch_model 
                    
                  }
                }
              }
            }
            
            # select optimal model according to AIC / BIC
              # initialize criterion value and index for best model
                aic_selected_model = c(10^10,0) #
                bic_selected_model = c(10^10,0)
                names(aic_selected_model)  = c("minimum_AIC","model_position")
                names(bic_selected_model)  = c("minimum_bic","model_position")
              
              # find model with lowest value for criterion
              for (i in 1:length(garch_model_selection)) {
                 if (aic_selected_model[1] >garch_model_selection[[i]]$model_evaluation$aic_model){
                    aic_selected_model[1:2] = c(garch_model_selection[[i]]$model_evaluation$aic_model, i)
                 }
                 if (bic_selected_model[1] >garch_model_selection[[i]]$model_evaluation$bic_model){
                    bic_selected_model[1:2] = c(garch_model_selection[[i]]$model_evaluation$bic_model, i)
                 }
              }
              
              # select model with lowest value for criterion
              print("Selected Models")
              if (aic_selected_model[2] == bic_selected_model[2]){
                print("Same Model selected by AIC and BIC")
                selected_model = garch_model_selection[[aic_selected_model[2]]]
                print("Model Specification")
                print(selected_model$model_specification)
              }
              if (aic_selected_model[2] != bic_selected_model[2]){
                print("Different Model selected by AIC and BIC")
                print("Select model manually or via different criterion")
                #selected_model = garch_model_selection[[aic_selected_model[2]]]
                #selected_model = garch_model_selection[[bic_selected_model[2]]]
                #print("Model Specification")
                #print(selected_model$model_specification)
              }
              
              all_selected_model[[data_iter]] = selected_model
          }
            
          # save estimated GARCH-model
          saveRDS(all_selected_model, file = paste0(outpathModels,"univariate_garchs_full_sample.rds"))
          
          
          ####Univariate Garch Opt function ####
      
          # ug_spec= ugarchspec(variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(0, ), include.mean = TRUE), distribution.model="norm")
          # ug_fit= ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
          # ug_fit

          opt_garch<- function(ar,ma,returns) {
            AIC <- matrix(data=NA,nrow=2,ncol=2)
            for (j in 1:3) {
              for (i in 1:3) {
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
          
