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
  
  garch_data_ts_r <- read_csv("data/data_outliers_1_with_values.csv",
                                          col_types = cols(DGS3MO = col_double(),
                                                           MOEX = col_double(), RUBEUR_val = col_double(),
                                                           SPX = col_double(), SPX_val = col_double(),
                                                           yc_1 = col_double(), yc_10 = col_double(),
                                                           yc_15 = col_double(), yc_2 = col_double(),
                                                           yc_20 = col_double(), yc_3 = col_double(),
                                                           yc_30 = col_double(), yc_5 = col_double(),
                                                           yc_7 = col_double()))
  cols = colnames(garch_data_ts_r)
  garch_data_ts_r = garch_data_ts_r[,2:ncol(garch_data_ts_r)]
  colnames(garch_data_ts_r) = cols[1:length(garch_data_ts_r)]

  print("Check match of colnames to values in dataset")
  print(head(garch_data_ts_r))
  
  # cut timeframe to after structural break
  garch_data_ts_r = filter(garch_data_ts_r, date>="2008-01-01")

  
  # insert mean for nas (for splits this does not have an impact)
  garch_data_ts_r = garch_data_ts_r %>%
    mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  garch_data_ts_r$date = as.Date(garch_data_ts_r$date)
  
 #garch_data_ts_r  = readRDS("output/univariateDescriptives/garch_data_ts_r.rds") 
 garch_data_ts_r_errors = as.xts(garch_data_ts_r[,c("rub_errors", "oil_errors")], order.by = garch_data_ts_r$date)
     
 
 
    
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
           returns_list_with_date = list(data.frame(date=index(garch_data_ts_r_errors), coredata(garch_data_ts_r_errors[,c("rub_errors")])),
                       data.frame(date=index(garch_data_ts_r_errors), coredata(garch_data_ts_r_errors[,c("oil_errors")]))) # save as data frame
           
           names(returns_list) = c("rub_all","oil_all") #!! rename if multiple timeframes are estimated
           names(returns_list_with_date) = c("rub_all","oil_all")

        # initialize selected model list for all univeriate series
          all_selected_model = vector("list", length = length(returns_list))
          names(all_selected_model) = names(returns_list)
          
        # estimate model for all univariate inputs
          for(data_iter in 1:length(returns_list)) {
            returns = returns_list[[data_iter]]
            returns_with_date = returns_list_with_date[[data_iter]]
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
                        garch_model = list(names(returns_list)[data_iter], returns, garch_coefs, model_specification, model_evaluation,returns_with_date)
                        names(garch_model) = c("series_name", "return_data", "garch_coefs", "model_specification", "model_evaluation","returns_with_date")
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
        
        base_split_variables = as.xts(cbind(epsilon, epsilon_sq)) # return and var as dataset split variables
        
        external_split_vars = c("rub_val","oil_val","MOEX_val","MOEX","SPX_val","SPX","RUBEUR_val","RUBEUR",
                                "yc_0.25","yc_1","yc_5","yc_10","DGS3MO","DGS10")
        additional_split_variables = as.xts(garch_data_ts_r[,external_split_vars], order.by = garch_data_ts_r$date)
        
        
        # get 2 lags of each main variable, 1 for others.
        max_lags = 3 # ! keep at 1+lags used. number lags for loop such that for all  split vars the same dataset is used (NA in the first observations otherwise)
        lag1 = lag(base_split_variables,1)
        lag2 = lag(base_split_variables,2)
        colnames(lag1) = paste0(colnames(base_split_variables),"_lag1")
        colnames(lag2) = paste0(colnames(base_split_variables),"_lag2")
        
        lag1_add = lag(additional_split_variables,1)
        colnames(lag1_add) = paste0(colnames(additional_split_variables),"_lag1")
        
        split_variables = as.data.frame(cbind(lag1 , lag2, lag1_add)) # data frame instead of time series
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
        
        # 4) Estimate and save GARCH for every terminal leaf (here only for RUBUSD)
        # input data
        # estimate model for each series and timeframe given
        number_timeframes = 1
        returns_list = list()
        returns_list_with_date = list()
        for (i in 1:length(subsamples_tree)) {
          returns_list[[i]]=subsamples_tree[[i]]$return # garch_data_ts_r only contains 1 dataframe. if multiple, then get all series for all timeframes in return_list
          returns_list_with_date[[i]]=subsamples_tree[[i]]
         }
        
        names(returns_list) = paste0("rub_subsample", seq(1:length(returns_list))) #!! rename if multiple timeframes are estimated
        names(returns_list_with_date) = paste0("rub_subsample", seq(1:length(returns_list_with_date)))
        
        # initialize selected model list for all submodels
        all_selected_model_tree = vector("list", length = length(returns_list))
        names(all_selected_model_tree) = names(returns_list)
        
        # estimate model for all univariate inputs
        for(data_iter in 1:length(returns_list)) {
          returns = returns_list[[data_iter]]
          returns_with_date = returns_list_with_date[[data_iter]]
          
          
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
          garch_model = list(names(returns_list)[data_iter], returns, garch_coefs, model_specification_tree, model_evaluation, returns_with_date,treeGarchResult$split_order_pruned,treeGarchResult$split_order)
          names(garch_model) = c("series_name", "return_data", "garch_coefs", "model_specification", "model_evaluation","returns_with_date","treeGarchResult", "treeGarchResultNoPruning")
          all_selected_model_tree[[data_iter]] = garch_model 
        }
        # save only tree
        saveRDS(all_selected_model_tree, file = paste0(outpathModels,"univariate_garchs_tree_rub.rds"))
        
        # investigate results
        for(i in 1:length(all_selected_model_tree)) {
          print(all_selected_model_tree[[i]]$series_name)
          print(all_selected_model_tree[[i]]$model_specification)
          print(all_selected_model_tree[[i]]$garch_coefs)
          print(all_selected_model_tree[[i]]$model_evaluation)
        }
        
        
      # OPTIONAL Model Saving for convinience: add simple model of other series split in same subsamples and save ----
        # ! only do if estimated on same time basis
        # ! need to adjust to number of splits done
        
        # tree and non-tree model - for one series
          # create list
          univ_rub_models_compare = list()
          univ_rub_models_compare$rub_m1_subsample1 = all_selected_model_tree$rub_subsample1
          univ_rub_models_compare$rub_m1_subsample2= all_selected_model_tree$rub_subsample2
          univ_rub_models_compare$rub_m2_subsample1 = all_selected_model$rub_all
          univ_rub_models_compare$rub_m2_subsample2 = all_selected_model$rub_all
          
          # change sample for non-tree model to tree subsamples
          univ_rub_models_compare$rub_m2_subsample1$return_data = univ_rub_models_compare$rub_m1_subsample1$return_data
          univ_rub_models_compare$rub_m2_subsample1$return_data_with_date = univ_rub_models_compare$rub_m1_subsample1$return_with_date 
          univ_rub_models_compare$rub_m2_subsample2$return_data = univ_rub_models_compare$rub_m1_subsample2$return_data
          univ_rub_models_compare$rub_m2_subsample2$return_data_with_date = univ_rub_models_compare$rub_m1_subsample2$return_with_date 
          
          # save list
          saveRDS(univ_rub_models_compare, file = paste0(outpathModels,"univ_rub_models_for_comparision.rds"))
          
        # tree and non-tree model - for different series
            # first match dates to get return series
            match_sample1 = select(inner_join(all_selected_model_tree$rub_subsample1$returns_with_date, all_selected_model$oil_all$returns_with_date, by="date"), return, date)
            match_sample2 = select(inner_join(all_selected_model_tree$rub_subsample2$returns_with_date, all_selected_model$oil_all$returns_with_date, by="date"), return, date)
            # match_sample3 = select(inner_join(all_selected_model_tree$rub_subsample3$returns_with_date, all_selected_model$oil_all$returns_with_date, by="date"), return, date)
            # match_sample4 = select(inner_join(all_selected_model_tree$rub_subsample4$returns_with_date, all_selected_model$oil_all$returns_with_date, by="date"), return, date)

            # create list
              univ_oil_rub_models_combined = list()
              univ_oil_rub_models_combined$rub_subsample1 = all_selected_model_tree$rub_subsample1
              univ_oil_rub_models_combined$rub_subsample2= all_selected_model_tree$rub_subsample2
              univ_oil_rub_models_combined$oil_subsample1 = all_selected_model$oil_all
              univ_oil_rub_models_combined$oil_subsample2 = all_selected_model$oil_all

            # change sample for non-tree model to tree subsamples
              univ_oil_rub_models_combined$rub_subsample1$return_data = select(match_sample1, -date)
              univ_oil_rub_models_combined$rub_subsample1$return_data_with_date = match_sample1
              univ_oil_rub_models_combined$rub_subsample2$return_data = select(match_sample2, -date)
              univ_oil_rub_models_combined$rub_subsample2$return_data_with_date = match_sample2
              
              univ_oil_rub_models_combined$oil_subsample1$return_data = select(match_sample1, -date)
              univ_oil_rub_models_combined$oil_subsample1$return_data_with_date = match_sample1
              univ_oil_rub_models_combined$oil_subsample2$return_data = select(match_sample2, -date)
              univ_oil_rub_models_combined$oil_subsample2$return_data_with_date = match_sample2
              
              # save list
              saveRDS(univ_oil_rub_models_combined, file = paste0(outpathModels,"univ_oil_rub_tree_models_combined_for_dcc.rds"))
 