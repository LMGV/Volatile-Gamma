# univariate garch models

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
  library(lubridate) 
  library(forecast) 
  library(tseries)
  library(xts)
  library(quantmod)
  library(fGarch)
  library(rugarch)
  library(rmgarch)
  library(psych)
  library(MASS)
  
  filter <- dplyr::filter
  select <- dplyr::select
  setwd("~/GitHub/Volatile-Gamma") # setwd
  
  source("scripts/functions.R") # functions
  source("scripts/garchFunction.R") # functions
  outpathDescriptive = "output/univariateDescriptives/"
  outpathModels =  "output/univariateModels/"
  
# Import Data ----
 garch_data_ts_r  = readRDS("output/univariateDescriptives/garch_data_ts_r.rds") # selected garch data after struc break analysis
 garch_data_ts_r_errors = garch_data_ts_r[,c("rub_errors", "oil_errors")]
     
 
   
       
       
       
       
     # forecasting
       # get model specification and parms:
       model_coefs= all_selected_model_tree$rub_subsample1$garch_coefs
       model_specif= all_selected_model_tree$rub_subsample1$model_specification
       returns_test = all_selected_model_tree$rub_subsample1$return_data[1:10]
       
       
       sigmasq = rep(var(returns_test),10)[1:max_lags_model] 
     # input: for each day one return vector + list model coefs + list model spec  
     # output: scalar variance forecast
       
       mu_coef = model_coefs[1]
       constant_coef = model_coefs[2]
       ma_coef = model_coefs[(1:model_specif$number_ma)+2]
       ar_coef = model_coefs[((model_specif$number_ma+1):(model_specif$number_ma+model_specif$number_ar))+2]
       th_coef = ifelse(model_specif$threshhold_included==T, model_coefs[(model_specif$number_ar+model_specif$number_ma+1)+2],0) #set threshhold to zero no TGARCH
        th_value = model_specif$th_value
      
       
     # input data frame
       max_lags_model = 3 # parameter that determines how many datapoints are loaded to forecast function
       
   
       
      
       

       # Note: timeindex: max_lags_model is t, (max_lags_model-1) is t-1 ...
       returns = returns_test[1:max_lags_model]
       
     # demeaned square return as variance proxy  
       epsilon  = returns -mu_coef$mu 
       epsilon_sq = epsilon^2
       
     # calc ma/ar/th parts of variance forecast
       ma_part = 0
         for (k in 1:length(ma_coef)) {
           ma_part = ma_part + ma_coef[k]*epsilon_sq[length(epsilon)+1-k]
         }
       ar_part = 0
         for (j in 1:length(ar_coef)) {
           ar_part = ar_part + ar_coef[j]*sigmasq[length(epsilon)+1-j]  # calc AR part
         }
       threshhold_part = th_coef*epsilon_sq[max_lags_model]* as.numeric(epsilon_sq[length(epsilon)]<=th_value)
       
     # variance estimation
       var_estim = as.numeric(constant_coef$omega + ar_part+ ma_part +threshhold_part)
       var_estim
       
       
       
            
          
      # Optimal fullsample GARCH model ----
         # possible model specifications
         ma_choices = seq(1:3)
         ar_choices = seq(1:3)
         threshhold_choices = c(T,F)
         th_value  = 0 # not optimized within fct
         data_threshhold = 0 # not implemented 
         distribution_choices =c("normal","t")

        # input data
         # estimate model for each series and timeframe given
           number_timeframes = 1
           returns_list=list(garch_data_ts_r_errors$rub_errors, garch_data_ts_r_errors$oil_errors) # garch_data_ts_r only contains 1 dataframe. if multiple, then get all series for all timeframes in return_list
               
           names(returns_list) = c("rub_all","oil_all") #!! rename if multiple timeframes are estimated

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
                        start_parms = c(0,0.1,  rep(0.1/ma,ma), rep(0.9/ar,ar)) # initialize parms. 
                        if(threshhold==T){
                          start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
                        }
                        if(distribution=="t"){
                          start_parms=  c(start_parms, 6) # keep df_t > 2 due to likelihood fct
                        }
                        
                        # get list for input specification
                          number_parms_estimated = length(start_parms)
                          model_specification = list(number_parms_estimated,ma,ar,threshhold, distribution,th_value, data_threshhold,start_parms)
                          names(model_specification) = c("number_parms_estimated","number_ma","number_ar","threshhold_included", "distribution","th_value","data_threshhold","start_parms")
                          
                        # estimate GARCH model for given specification (minimize negative loglikelihood)
                          opt_parms= nlm(garchEstimation,model_specification$start_parms,
                                         returns = returns,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                                         threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = model_specification$data_threshhold,
                                         distribution=model_specification$distribution,
                                         print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
                        
                        # get model parameters
                          # get same sames as in DCC function
                          names_parms =  c("mu", "omega", paste0("alpha",seq(1:ma)), paste0("beta",seq(1:ar)))
                          if(model_specification$threshhold==T){
                            names_parms=  c(names_parms, "eta11") # set asymmetry parameter to 0 
                          }
                          if(model_specification$distribution=="t"){
                            names_parms=  c(names_parms, "shape") # keep df_t > 2 due to likelihood fct
                          }
                          
                          garch_coefs = as.data.frame(t(c(opt_parms$estimate[1], opt_parms$estimate[2:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma):length(opt_parms$estimate)])))
                          colnames(garch_coefs) = names_parms
                        
                        # model evaluation
                            # stationarity
                            sum_coefs = sum(garch_coefs[,3:(2+ar+ma)])
                            if(model_specification$threshhold==T){
                              sum_coefs=  sum_coefs + sum(returns<=th_value)/(length(returns))*garch_coefs$eta11 # adjust if threshhold is active
                            }
                            
                            # model selection
                            loglik_model =-opt_parms$minimum
                            aic_model = my_aic(loglik_model, model_specification$number_parms_estimated)
                            bic_model = my_bic(loglik_model, model_specification$number_parms_estimated, length(returns))
                            
                            # save model evaluation in list  
                            model_evaluation = list(sum_coefs,loglik_model, aic_model, bic_model)
                            names(model_evaluation) = c("sum_ar_ma_coefs","log_lik","aic_model","bic_model")
                        
                       # add model to model selection list
                        garch_model = list(names(returns_list)[data_iter], returns, garch_coefs, model_specification, model_evaluation)
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
              if (aic_selected_model[2] == bic_selected_model[2]){
                print("Same Model selected by AIC and BIC")
                selected_model = garch_model_selection[[aic_selected_model[2]]]
                print("Model Specification")
                print(selected_model$model_specification)
              }
              if (aic_selected_model[2] != bic_selected_model[2]){
                print("Different Model selected by AIC and BIC")
                print("Select Model According to BIC")
                selected_model = garch_model_selection[[bic_selected_model[2]]]
                #selected_model = garch_model_selection[[bic_selected_model[2]]]
                #print("Model Specification")
                #print(selected_model$model_specification)
              }
              
              all_selected_model[[data_iter]] = selected_model
          }
        # save estimated GARCH-model
        saveRDS(all_selected_model, file = paste0(outpathModels,"univariate_garchs_full_sample.rds"))
        
        # investigate results
        for(i in 1:length(all_selected_model)) {
          print(all_selected_model[[i]]$series_name)
          print(all_selected_model[[i]]$model_specification)
          print(all_selected_model[[i]]$garch_coefs)
          print(all_selected_model[[i]]$model_evaluation)
        }

        
# Tree GARCH (1,1) ----
  # 1) define data inputs (settings: only for RUBUSD, can easily be extended)
    # define possible split variables
      # past returns, epsilons, variances of own process and other processes
        means = colMeans(garch_data_ts_r_errors)
        epsilon = sweep(garch_data_ts_r_errors,2,means)
        epsilon_sq = epsilon^2
        colnames(epsilon) = paste0(colnames(garch_data_ts_r_errors),"_epsilon")
        colnames(epsilon_sq) = paste0(colnames(garch_data_ts_r_errors),"_epsilon_sq")
        
        base_split_variables = as.xts(cbind(epsilon, epsilon_sq)) #dataset split variables
        
        # get 2 lags of each variable. In VAR tests, dependencies were not consistent above the 2nd lag. Are excluded to get parsimonious computationally feasible model
        max_lags = 3 # ! keep at 1+lags used. number lags for loop such that for all  split vars the same dataset is used (NA in the first observations otherwise)
        lag1 = lag(base_split_variables,1)
        lag2 = lag(base_split_variables,2)
        colnames(lag1) = paste0(colnames(base_split_variables),"_lag1")
        colnames(lag2) = paste0(colnames(base_split_variables),"_lag2")
        
        split_variables = as.data.frame(cbind(lag1 , lag2)) # data frame instead of time series
        split_variables$date = rownames(split_variables)
        
        vector_quantiles = seq(1, 7)*0.125 # quantiles as threshholds
        
        # list for splitting variables
        list_split_variables = colnames(split_variables)[colnames(split_variables) !="date"] # cols as splitting variables
        
        # return series used
        returns = as.data.frame(garch_data_ts_r_errors$rub_errors)
        returns$date = rownames(returns)
        colnames(returns) = c("return","date")
        
        # 2) define model specifications for tree - GARCH (1,1)
        ma = 1
        ar = 1
        threshhold = F
        th_value  = 0 # not optimized within fct
        data_threshhold = 0 # not implemented 
        distribution =c("t") # take normal. convergence issues with t in small subsamples of tree  
        start_parms = c(0,0.1,  rep(0.1/ma,ma), rep(0.9/ar,ar)) # initialize parms. 
        if(threshhold==T){
          start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
        }
        if(distribution=="t"){
          start_parms=  c(start_parms, 10) # keep df_t > 2 due to likelihood fct
        }
        number_parms_estimated = length(start_parms)
        model_specification_tree = list(number_parms_estimated,ma,ar,threshhold, distribution,th_value, data_threshhold,start_parms)
        names(model_specification_tree) = c("number_parms_estimated","number_ma","number_ar","threshhold_included", "distribution","th_value","data_threshhold","start_parms")
        
        
        
        # 2) build tree
        treeGarchResult =   buildAndPruneTree(returns, split_variables, list_split_variables, model_specification_tree, max_lags)  
        
        
        # 3) get subsamples tree
        subsamples_tree  = collectSubsamplesTree(returns, split_variables, treeGarchResult$split_order_pruned)
        
        # 4) optimal GARCH for every terminal leaf (here only for RUBUSD)
        # input data
        # estimate model for each series and timeframe given
        number_timeframes = 1
        returns_list = list()
        for (i in 1:length(subsamples_tree)) {
          returns_list[[i]]=subsamples_tree[[i]]$return # garch_data_ts_r only contains 1 dataframe. if multiple, then get all series for all timeframes in return_list
        }
        
        names(returns_list) = paste0("rub_subsample", seq(1:length(returns_list))) #!! rename if multiple timeframes are estimated
        
        # initialize selected model list for all submodels
        all_selected_model_tree = vector("list", length = length(returns_list))
        names(all_selected_model_tree) = names(returns_list)
        
        # estimate model for all univariate inputs
        for(data_iter in 1:length(returns_list)) {
          returns = returns_list[[data_iter]]
          
          
          # estimate GARCH model for given specification_tree (minimize negative loglikelihood)
          opt_parms= nlm(garchEstimation,model_specification_tree$start_parms,
                         returns = returns,  ma = model_specification_tree$number_ma, ar = model_specification_tree$number_ar, 
                         threshhold = model_specification_tree$threshhold_included, th_value = model_specification_tree$th_value, data_threshhold = model_specification_tree$data_threshhold,
                         distribution=model_specification_tree$distribution,
                         print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
          
          # get model parameters
          # get same names as in DCC function
          names_parms =  c("mu", "omega", paste0("alpha",seq(1:ma)), paste0("beta",seq(1:ar)))
          if(model_specification_tree$threshhold==T){
            names_parms=  c(names_parms, "eta11") # set asymmetry parameter to 0 
          }
          if(model_specification_tree$distribution=="t"){
            names_parms=  c(names_parms, "shape") # keep df_t > 2 due to likelihood fct
          }
          
          garch_coefs = as.data.frame(t(c(opt_parms$estimate[1], opt_parms$estimate[2:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma):length(opt_parms$estimate)])))
          colnames(garch_coefs) = names_parms
          
          # model evaluation
          # stationarity
          sum_coefs = sum(garch_coefs[,3:(2+ar+ma)])
          if(model_specification_tree$threshhold==T){
            sum_coefs=  sum_coefs + sum(returns<=th_value)/(length(returns))*garch_coefs$eta11 # adjust if threshhold is active
          }
          
          # model selection
          loglik_model =-opt_parms$minimum
          aic_model = my_aic(loglik_model, model_specification_tree$number_parms_estimated)
          bic_model = my_bic(loglik_model, model_specification_tree$number_parms_estimated, length(returns))
          
          # save model evaluation in list  
          model_evaluation = list(sum_coefs,loglik_model, aic_model, bic_model)
          names(model_evaluation) = c("sum_ar_ma_coefs","log_lik","aic_model","bic_model")
          
          # add model to model selection list
          garch_model = list(names(returns_list)[data_iter], returns, garch_coefs, model_specification_tree, model_evaluation, treeGarchResult$split_order_pruned)
          names(garch_model) = c("series_name", "return_data", "garch_coefs", "model_specification", "model_evaluation","treeGarchResult")
          all_selected_model_tree[[data_iter]] = garch_model 
        }
        # save only tree
        saveRDS(all_selected_model_tree, file = paste0(outpathModels,"univariate_garchs_tree_rub.rds"))
        # add optimal 
        
        # investigate results
        for(i in 1:length(all_selected_model)) {
          print(all_selected_model[[i]]$series_name)
          print(all_selected_model[[i]]$model_specification)
          print(all_selected_model[[i]]$garch_coefs)
          print(all_selected_model[[i]]$model_evaluation)
        }
        
        
        
# REMOVE LATER ----
      
        
        acf(returns_list$rub_timeframe3$rub_errors^2)
        pacf(returns_list$rub_timeframe3$rub_errors^2)
        test_data = returns_list$oil_timeframe3$oil_errors
        algorithm = c("nlminb", "lbfgsb", "nlminb+nm", "lbfgsb+nm")
        garchFit(~ garch(1, 1), data = test_data, trace = F, cond.dist =  "std")
        
        
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
          
          ####Uni-Garch Specs Oil and RUB#
          #output_oil <- opt_garch(0,0, test_data)
          output_rub <- opt_garch(0,0, test_data)
          
          #Garch summary Oil
          output_oil[["fit"]]
          
          #Garch summary Oil
          output_rub[["fit"]]
          
          ####Bivar Garch#
          # cor(ts_r$rub_errors,ts_r$oil_errors)
          # var(ts_r)
          # var(ts_r^2)
          # cor(ts_r$rub_errors^2,ts_r$oil_errors^2)
          
