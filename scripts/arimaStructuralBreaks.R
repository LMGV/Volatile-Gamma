setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub")
library(lubridate) 
source("scripts/functions.R") # functions
source("scripts/functionsARIMA.R")
outpathDescriptive = "output/structARIMA/" 


data.r <- read.table("data/data_outliers_1_with_values.csv", sep = ",")
data.r$date <- as.Date(data$date, "%Y-%m-%d")
loglik.full <- read.table(file = "data/ARIMA_LL.txt") #LL values of the chosen ARIMA models for oil and rub, respectively
estimated.p <- read.table(file = "data/ARIMA_p.txt")
estimated.q <- read.table(file = "data/ARIMA_q.txt")


number_restrictions.oil = estimated.p[1, 1] + estimated.q[1, 1] + 1 # coefficients of AR and MA parts, plus a constant
number_restrictions.rub = estimated.p[2, 1] + estimated.q[2, 1] + 1

# test for structural breaks ----
  # model specification for testing struc breaks
    
#@mila: dont need this if you choose optimal model within the function
    
  # test for 2 complete models vs 1 complete model
    # @mila number_restrictions = parms unrestricted model - parms restricted model.
    # @Mila retricted = same model for entire timeframe, unrestric = 2 models (for before and after). Probably just estimate complete new model instead of only certain parms (easier to implement)
   
    # @Mila: if you
    significance_level = 0.01

    # input data
    returns_list=list(data.r$oil, data.r$rub)
    names(returns_list) = c("oil","rub")
    
    # Choice univariate return series
        # test all series and checked that main struc breaks for all series are covered. 
        # our heuristic method: set same breaks for both series
        # 1. set ruble structural breaks (key series of interest) 
        # 2. check performance of those breaks for oil.
    
    # loop over structural break grid: 
    # procedure: iteratively pick statistically and economically plausible structural break 
    # check break each quarter
    min_time_estim = 1
    start_date = as.Date(as.yearqtr(returns$date[1], format="%Y-%m-%d"))
    end_date = as.Date(as.yearqtr(returns$date[nrow(returns)], format="%Y-%m-%d"))
    grid_struct_breaks = seq(start_date+years(min_time_estim), end_date-years(min_time_estim), by = "quarter") # grid of structural breaks. atleast min_time_estim years before and after structural break for estimation
    
    ### searching for structural breaks in oil series
    returns.oil = data.frame(date = data.r$date, ret = data.r$oil)
    table_struc_break1_arima.oil = find_structural_break_arima(returns.oil, loglik.full[1,1], grid_struct_breaks, number_restrictions.oil, 
                                                               significance_level)
    write.csv(table_struc_break1_arima.oil, file = paste0(outpathDescriptive,"table_struc_break1_arima_oil.csv"))
    
    ### searching for structural breaks in rub series
    returns.rub = data.frame(date = data.r$date, ret = data.r$rub)
    table_struc_break1_arima.rub = find_structural_break_arima(returns.rub, loglik.full[2,1], grid_struct_breaks, number_restrictions.rub, 
                                                               significance_level)
    write.csv(table_struc_break1_arima.rub, file = paste0(outpathDescriptive,"table_struc_break1_arima_rub.csv"))
        
          
find_structural_break_arima = function(returns, LL_full, grid_struct_breaks, number_restrictions, significance_level) {
  # initialize test-result table
  struc_break_test_results = as.data.frame(matrix(nrow = length(grid_struct_breaks), ncol = 7,))
  colnames(struc_break_test_results) = c("break_date","test_result_word","test_result","p_value","test_stat" , "logLikUR", "logLikR")
  struc_break_test_results$break_date = grid_struct_breaks
  
  # estimate 2 Arima models for before and after structural breaks in grid 
  for (break_iter in 1:length(grid_struct_breaks)) {
    sample_before = returns[returns[[1]] <  grid_struct_breaks[break_iter],][[2]]
    sample_after  = returns[returns[[1]] >= grid_struct_breaks[break_iter],][[2]]
    
    # estimate Arima model before and after potential breakpoint
    model_before = run.ARIMA.fit(sample_before)
    model_after  = run.ARIMA.fit(sample_after)
    
    # number_restrictions
    
    # LR-test
    struc_break_test_results[break_iter,2:7] = supLikelihoodTest(model_before$LogLikelihoodvalue + model_after$LogLikelihoodvalue,
                                                                 LL_full, number_restrictions, significance_level) 
  }
  return(struc_break_test_results)
}  

