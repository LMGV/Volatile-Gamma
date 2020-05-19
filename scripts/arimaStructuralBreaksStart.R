source("scripts/functions.R") # functions
outpathDescriptive = "output/SOMETHING" ##### specify output path for tables


# test for structural breaks ----
  # model specification for testing struc breaks
    
#@mila: dont need this if you choose optimal model within the function
    
  # test for 2 complete models vs 1 complete model
    # @mila number_restrictions = parms unrestricted model - parms restricted model.
    # @Mila retricted = same model for entire timeframe, unrestric = 2 models (for before and after). Probably just estimate complete new model instead of only certain parms (easier to implement)
   
    # @Mila: if you
    significance_level = 0.01

    # input data
    returns_list=list(ts_r$rub_errors, ts_r$oil_errors)
    names(returns_list) = c("rub","oil")
    
    # Choice univariate return series
        # test all series and checked that main struc breaks for all series are covered. 
        # our heuristic method: set same breaks for both series
        # 1. set ruble structural breaks (key series of interest) 
        # 2. check performance of those breaks for oil.
        for (series_iter in 1:length(returns_list))
        {
          # set data in loop
          returns = returns_list[[series_iter]] 
          
          # loop over structural break grid: 
          # procedure: iteratively pick statistically and economically plausible structural break 
          # check break each quarter
          min_time_estim = 1
          start_date = as.Date(as.yearqtr(time(returns)[1], format="%Y-%m-%d"))
          end_date = as.Date(as.yearqtr(time(returns)[length(returns)], format="%Y-%m-%d"))
          grid_struct_breaks = seq(start_date+years(min_time_estim), end_date-years(min_time_estim), by = "quarter") # grid of structural breaks. atleast min_time_estim years before and after structural break for estimation
          
          #break1
          # @Mila: This is the annoying part. Need to rewrite this function. see end of this file
          table_struc_break1_arima = find_structural_break(returns,grid_struct_breaks, start_parms, model_specification, number_restrictions,significance_level)
          write.csv(table_struc_break1_arima, file = paste0(outpathDescriptive,"table_struc_break1_arima_",names(returns),".csv"))
          
          # manually select by p_values and chart
          # @ Mila: look at the result table of the test and pick reasonable struc break
          # @ MILA: I used this break for garch. If its reasonable then also use it, otherwise specify on your own
          struc_break1 = as.Date("2008-01-01") 
          print("Set struc break for both series")
          print(struc_break1)
          
          # break2
          returns_break2_1  =returns[time(returns) < struc_break1]
          returns_break2_2 =returns[time(returns) >=struc_break1]
          
          grid_struct_breaks2_1 = grid_struct_breaks[(grid_struct_breaks>=start_date+years(min_time_estim)) & (grid_struct_breaks <= (struc_break1-years(min_time_estim)))]
          grid_struct_breaks2_2 = grid_struct_breaks[(grid_struct_breaks<=end_date-years(min_time_estim)) & (grid_struct_breaks >= (struc_break1+years(min_time_estim)))]
          
          # @Mila: change fct
          table_struc_break2_1_arima = find_structural_break(returns_break2_1,grid_struct_breaks2_1, start_parms, model_specification, number_restrictions, significance_level)
          table_struc_break2_2_arima = find_structural_break(returns_break2_2,grid_struct_breaks2_2, start_parms, model_specification, number_restrictions,significance_level)
          
          write.csv(table_struc_break2_1_arima, file = paste0(outpathDescriptive,"table_struc_break2_1_arima_",names(returns),".csv"))
          write.csv(table_struc_break2_2_arima, file = paste0(outpathDescriptive,"table_struc_break2_2_arima_",names(returns),".csv"))
          
          
          # manually select by looking at table
          # @ Mila: IN Garch we in the end only use the 2005 struc break and exclude the part before. 
          struc_break2_1 = as.Date('2005-01-01') # minimal p-value for rub. No breaks for oil, but will still apply it to both
          struc_break2_2 = as.Date('2014-04-01') # no rejection of null hypothesis for rub. rejection on '2014-04-01' for oil, p-value of rub at 0.08. apply break
          
          # break3
          returns_break3_1  =returns_break2_1[time(returns_break2_1) < struc_break2_1]
          returns_break3_2 =returns_break2_1[time(returns_break2_1) >=struc_break2_1]
          returns_break3_3  =returns_break2_2[time(returns_break2_2) < struc_break2_2]
          returns_break3_4 =returns_break2_2[time(returns_break2_2) >=struc_break2_2]
          
          grid_struct_breaks3_1 = grid_struct_breaks[(grid_struct_breaks>=start_date+years(min_time_estim)) & (grid_struct_breaks <= (struc_break2_1-years(min_time_estim)))]
          grid_struct_breaks3_2 = grid_struct_breaks[(grid_struct_breaks>=struc_break2_1+years(min_time_estim)) & (grid_struct_breaks <= (struc_break1-years(min_time_estim)))]
          grid_struct_breaks3_3 = grid_struct_breaks[(grid_struct_breaks>=struc_break1+years(min_time_estim)) & (grid_struct_breaks <= (struc_break2_2-years(min_time_estim)))]
          grid_struct_breaks3_4 = grid_struct_breaks[(grid_struct_breaks>=struc_break2_2+years(min_time_estim)) & (grid_struct_breaks <= (end_date-years(min_time_estim)))]
          
          # @Mila: change fct
          table_struc_break3_1_arima = find_structural_break(returns_break3_1,grid_struct_breaks3_1, start_parms, model_specification, number_restrictions, significance_level)
          table_struc_break3_2_arima = find_structural_break(returns_break3_2,grid_struct_breaks3_2, start_parms, model_specification, number_restrictions,significance_level)
          table_struc_break3_3_arima = find_structural_break(returns_break3_3,grid_struct_breaks3_3, start_parms, model_specification, number_restrictions, significance_level)
          table_struc_break3_4_arima = find_structural_break(returns_break3_4,grid_struct_breaks3_4, start_parms, model_specification, number_restrictions,significance_level)
          
          write.csv(table_struc_break3_1_arima, file = paste0(outpathDescriptive,"table_struc_break3_1_arima_",names(returns),".csv"))
          write.csv(table_struc_break3_2_arima, file = paste0(outpathDescriptive,"table_struc_break3_2_arima",names(returns),".csv"))
          write.csv(table_struc_break3_3_arima, file = paste0(outpathDescriptive,"table_struc_break3_3_arima",names(returns),".csv"))
          write.csv(table_struc_break3_4_arima, file = paste0(outpathDescriptive,"table_struc_break3_4_arima",names(returns),".csv"))
        }
        
        # save return for each timeperiod
        # only use struc_break1 and struc_break2_1 (struc_break2_2 is not sign for both. requires)
          # @mila: cut timeframes for different arima models according to structural breaks 
          # @ Mila: check for stationarity again

          ts_r_time1 = ts_r[(time(ts_r) < struc_break2_1)]
          ts_r_time2 = ts_r[(time(ts_r) >= struc_break2_1) & (time(ts_r)< struc_break1)]
          ts_r_time3 = ts_r[(time(ts_r) >= struc_break1)]
          
          ts_r_struc_break = list(ts_r_time1,ts_r_time2,ts_r_time3)
          names(ts_r_struc_break) = c("timeframe1", "timeframe2", "timeframe3")
          
          # @ Mila: !!after getting the residuals merge it together in one timeframe again
          
          

# @ Mila Rewrite this thing. give different name to fct
  # @ Mila  inputs to fct you need to change/redefine: start parms(you probably dont even need that?), model_specification(ommit if you choose opt model in loop), ...
          # @ Mila ...  number_restrictions (calculate within loop if you do auto model selection)
          
find_structural_break = function(returns,grid_struct_breaks, start_parms, model_specification, number_restrictions,significance_level) {
  # initialize test-result table
  struc_break_test_results = as.data.frame(matrix(nrow = length(grid_struct_breaks), ncol = 7,))
  colnames(struc_break_test_results) = c("break_date","test_result_word","test_result","p_value","test_stat" , "logLikUR", "logLikR")
  struc_break_test_results$break_date = grid_struct_breaks
  
  # estimate Arima model for entire timeframe
     # @ Mila required output: log-likelihoodof arima full timeframe. If its negative log-lik then change sign in LR-test
  
  # estimate 2 Arima models for before and after structural breaks in grid 
  for (break_iter in 1:length(grid_struct_breaks)) {
    sample_before = returns[index(returns) < grid_struct_breaks[break_iter]]
    sample_after   = returns[index(returns) >= grid_struct_breaks[break_iter]]
    
    # estimate Arima model before and after potential breakpoint
    # @ Mila required output: log-likelihoodof arima before and arima after struc break. If its negative log-lik then change sign in LR-test
    
    # number_restrictions
    # @ Mila: calculate number_restrictions = parms_unrestricted - parms_restricted
    
    # LR-test
    struc_break_test_results[break_iter,2:7] = supLikelihoodTest(-opt_parms_before$minimum - opt_parms_after$minimum,-opt_parms$minimum, number_restrictions,significance_level) 
  }
  return(struc_break_test_results)
}  

