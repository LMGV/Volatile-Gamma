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
     
 
 # Tree GARCH (1,1) ----
            returns = as.data.frame(garch_data_ts_r_errors$rub_errors)
            returns$date = rownames(returns)
            colnames(returns) = c("return","date")

              # define possible split variables
                # past returns, epsilons, variances of own process and other processes
                   means = colMeans(garch_data_ts_r_errors)
                   epsilon = sweep(garch_data_ts_r_errors,2,means)
                   epsilon_sq = epsilon^2
                   colnames(epsilon) = paste0(colnames(garch_data_ts_r_errors),"_epsilon")
                   colnames(epsilon_sq) = paste0(colnames(garch_data_ts_r_errors),"_epsilon_sq")
                   
                   base_split_variables = as.xts(cbind(epsilon, epsilon_sq)) #dataset split variables
                 
                # get 2 lags of each variable. In VAR tests, dependencies were not consistent above the 2nd lag. Are excluded to get parsimonious computationally feasible model
                   max_lags = 2 # number lags for loop such that for all  split vars the same dataset is used (NA in the first observations otherwise)
                   lag1 = lag(base_split_variables,1)
                   lag2 = lag(base_split_variables,2)
                   colnames(lag1) = paste0(colnames(base_split_variables),"_lag1")
                   colnames(lag2) = paste0(colnames(base_split_variables),"_lag2")
                   
                   split_variables = as.data.frame(cbind(lag1 , lag2)) # data frame instead of time series
                   split_variables$date = rownames(split_variables)

                   vector_quantiles = seq(1, 7)*0.125 # quantiles as threshholds
                 
                 # list for splitting variables
                   list_split_variables = colnames(split_variables)[colnames(split_variables) !="date"] # cols as splitting variables
                 
             # dataframes to collect subsamples and split information
                   
               leaf_history = list() # leaf history contains for each split level: the sample, the "selected_split", the history
                   
               history_names = c("parent","own") # each leaf gets a history vector with parent number and a own number for clear assignment
               
               active_leafs = list() # collect active samples which can be used to split

               history_start = c(0,0) # starting sample gets itself as parent. Stop 
               leaf_history[[1]] = list(returns, NULL, history_start)
               
               names_sample = c("sample","parent_split")
               
               # first split
                  # run split function for both nodes and choose best split according to likelihood improvement 
                   split_1 = find_split(returns, split_variables, list_split_variables, max_lags)
                  # execute split (since only one choosen split)
                    data_sample1 = list(split_1$subsample_lower, split_1)
                    data_sample2 = list(split_1$subsample_higher,split_1)
                    names(data_sample1) = names_sample
                    names(data_sample2) = names_sample
                    active_leafs = list(data_sample1, data_sample2) # overwrite in first split
                  
                  # leaf_history[[2]] = list(split_1$subsample_lower, split_1$selected_split, c(0,1))
                  # leaf_history[[3]] = list(split_1$subsample_higher, split_1$selected_split, c(0,2))
                   
                # second split
                  # run split function for both nodes and choose best split according to likelihood improvement
                  split_2_1 = find_split(active_leafs[[1]]$sample, split_variables, list_split_variables, max_lags)
                  split_2_2 =  find_split(active_leafs[[2]]$sample, split_variables, list_split_variables, max_lags)
                  
                  # split if best likelihood improvement and update active leafs
                    if(split_2_1$selected_split$improvementLogLik >= split_2_2$selected_split$improvementLogLik) {
                      split_2 = split_2_1
                      
                      data_sample3 = list(split_2$subsample_lower, split_2)
                      data_sample4 = list(split_2$subsample_higher,split_2)
                      data_sample5 = list(split_1$subsample_higher,split_1)
                      names(data_sample1) = names_sample
                      names(data_sample2) = names_sample
                      
                      active_leafs[[1]] = active_leafs[[2]] 
                      # leaf_history[[4]] = list(split_2$subsample_lower, split_2$selected_split, c(1,3))
                      # leaf_history[[5]] = list(split_2$subsample_higher, split_2$selected_split, c(1,4))
                      # leaf_history[[6]] = list(split_1$subsample_higher, split_1$selected_split, c(0,2))
                      
                    } else {
                      split_2 = split_2_2
                      active_leafs[[1]] = active_leafs[[1]] 
                      # leaf_history[[4]] = list(split_2$subsample_lower, split_2$selected_split, c(2,3))
                      # leaf_history[[5]] = list(split_2$subsample_higher, split_2$selected_split, c(2,4))
                      # leaf_history[[6]] = list(split_1$subsample_lower, split_1$selected_split, c(0,1))
                    }
                  
                  # execute split
                    active_leafs[[2]] =  split_2$subsample_lower
                    active_leafs[[3]] =  split_2$subsample_higher

                # third split
                  # run split function for all nodes and choose best split according to likelihood improvement
                  split_3_1 = find_split(active_leafs[[1]], split_variables, list_split_variables, max_lags)
                  split_3_2 = find_split(active_leafs[[2]], split_variables, list_split_variables, max_lags)
                  split_3_3 = find_split(active_leafs[[3]], split_variables, list_split_variables, max_lags)
                  
                  # split if best likelihood improvement and update active leafs
                    if((split_3_1$selected_split$improvementLogLik >= split_3_2$selected_split$improvementLogLik) & (split_3_1$selected_split$improvementLogLik >= split_3_3$selected_split$improvementLogLik)) {
                      split_3 = split_3_1
                      active_leafs[[1]] = active_leafs[[2]]
                      active_leafs[[2]] = active_leafs[[3]]
                      
                      leaf_history[[7]] = list(split_2$subsample_lower, split_2$selected_split, c(1,3))
                      leaf_history[[8]] = list(split_2$subsample_higher, split_2$selected_split, c(1,4))
                      leaf_history[[9]] = list(split_1$subsample_lower, split_1$selected_split, c(3,1))
                      leaf_history[[10]] = list(split_1$subsample_lower, split_1$selected_split, c(4,1))
                      
                    } else if ((split_3_2$selected_split$improvementLogLik > split_3_1$selected_split$improvementLogLik) & (split_3_2$selected_split$improvementLogLik >= split_3_3$selected_split$improvementLogLik)){
                      split_3 = split_3_2
                      active_leafs[[1]] = active_leafs[[1]]
                      active_leafs[[2]] = active_leafs[[3]]
                    } else {
                      split_3 = split_3_3
                      active_leafs[[1]] = active_leafs[[1]]
                      active_leafs[[2]] = active_leafs[[2]]
                    }
                  # execute split
                  active_leafs[[3]] =  split_3$subsample_lower
                  active_leafs[[4]] =  split_3$subsample_higher
                  
                  # save split history
                  split_history= list(split_1, split_2, split_3)
                  names(split_history) = c("split_1", "split_2", "split_3")
                  
             
                print(split_history$split_1$selected_split)
                print(split_history$split_2$selected_split)
                print(split_history$split_3$selected_split)
                  

                 
              # step 3) prune: choose subtree that minimized sum AIC
                 
            

               
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
                          model_specification = list(number_parms_estimated,ma,ar,threshhold, distribution,th_value)
                          names(model_specification) = c("number_parms_estimated","number_ma","number_ar","threshhold_included", "distribution","th_value")
                          
                        # estimate GARCH model for given specification (minimize negative loglikelihood)
                          opt_parms= nlm(garchEstimation,start_parms,
                                         returns = returns,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                                         threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                                         distribution=model_specification$distribution,
                                         print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
                        
                        # get model parameters
                          # get same sames as in DCC function
                          names_parms =  c("mu", "omega", paste0("alpha",seq(1:ma)), paste0("beta",seq(1:ar)))
                          if(threshhold==T){
                            names_parms=  c(names_parms, "eta11") # set asymmetry parameter to 0 
                          }
                          if(distribution=="t"){
                            names_parms=  c(names_parms, "shape") # keep df_t > 2 due to likelihood fct
                          }
                          
                          garch_coefs = as.data.frame(t(c(opt_parms$estimate[1], opt_parms$estimate[2:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma):length(opt_parms$estimate)])))
                          colnames(garch_coefs) = names_parms
                        
                        # model evaluation
                            # stationarity
                            sum_coefs = sum(garch_coefs[,3:(2+ar+ma)])
                            if(threshhold==T){
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

        
          ####Univariate Garch Opt function ####
      
          # ug_spec= ugarchspec(variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(0, ), include.mean = TRUE), distribution.model="norm")
          # ug_fit= ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
          # ug_fit
        
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
          
