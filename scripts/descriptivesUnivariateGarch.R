# Univeriate Descriptives for Squares

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
library(FinTS)

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

# Purely descriptive volatility plot ----
  # plot 30d ma daily volatility over time
  # proxy: mean of 30 day epsilon^2 = (value-mean(30days))^2
    volatility_plot_data = as.data.frame(ts_r)
    
    volatility_plot_data = arrange(volatility_plot_data, index(volatility_plot_data)) %>%
      mutate(Date = time(ts_r),
             sd_rub = sqrt(252*rollapply((rub_errors-rollapply(rub_errors,360,mean, na.rm = TRUE, align = "right", fill = NA))^2, 30,mean, na.rm = TRUE, align = "right", fill = NA)),
             sd_oil = sqrt(252*rollapply((oil_errors-rollapply(oil_errors,360,mean, na.rm = TRUE, align = "right", fill = NA))^2, 30,mean, na.rm = TRUE, align = "right", fill = NA))) %>%
      filter(is.na(sd_rub)==F & is.na(sd_oil)==F)
    
    summary(volatility_plot_data)
    title = "RUBUSD_and_Oil_Volatility"
    xlab = "Time"
    ylab= "Annual Volatility (30d-MA of squared deviation)"
    y1 = volatility_plot_data$sd_rub
    y2 = volatility_plot_data$sd_oil
    x = volatility_plot_data$Date
    names_y = c("RUB/USD","Oil")
    line_plot_multiple(title, outpath = outpathDescriptive, x,xlab, ylab, names_y, y_percent=T, y_discrete=F, legend=T, y1, y2)


# Autocorrelation of squares ----
  # ACF / PACF r^2 
  ts_r2 <- ts_r^2
  acf(ts_r2$oil_errors, lag.max = 30, plot = TRUE)
  acf(ts_r2$rub_errors, lag.max = 30, plot = TRUE)  
  pacf(ts_r2$oil_errors, lag.max = 30, plot = TRUE)
  pacf(ts_r2$rub_errors, lag.max = 30, plot = TRUE)  

  # Lagrange Multiplier test for presence of ARCH
  ArchTest (ts_r2$oil_errors, lags=30, demean = TRUE)
  ArchTest (ts_r2$rub_errors, lags=30, demean = TRUE)
  # result: for entire series, ARCH effects are present
  
# Structural breaks ----
  # test for structural breaks ----
    # model specification for testing struc breaks
      # GARCH 1,1 with t-distrib (TGarch was only useful for oil, will probably yield similar structural breaks)
        ma = 1
        ar = 1
        threshhold = F
        th_value  = 0 # not optimized within fct
        data_threshhold = 0 # not implemented 
        distribution =c("t")
    
      # starting parms
        start_parms = c(0,0.1,  rep(0.1/ma,ma), rep(0.9/ar,ar)) # initialize parms. 
        if(threshhold==T){
          start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
        }
        if(distribution=="t"){
          start_parms=  c(start_parms, 4) # keep df_t > 2 due to likelihood fct
        }
    
    
      # combine model specification
        number_parms_estimated = length(start_parms)
        model_specification = list(number_parms_estimated,ma,ar,threshhold, distribution,th_value)
        names(model_specification) = c("number_parms_estimated","number_ma","number_ar","threshhold_included", "distribution","th_value")
    
# STRUCTURAL BREAKS ----
  # output: datasets for structural breaks
  # take simple GARCH 1,1 model with t-errors to determine structural breaks
    # test for 2 complete models vs 1 complete model
    number_restrictions = length(start_parms) # number of restrictions is number parms 1 model
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
        table_struc_break1 = find_structural_break(returns,grid_struct_breaks, start_parms, model_specification, number_restrictions,significance_level)
        write.csv(table_struc_break1, file = paste0(outpathDescriptive,"table_struc_break1_",names(returns),".csv"))
      
        # manually select by p_values and chart
          struc_break1 = as.Date("2008-01-01") 
          print("Set struc break for both series")
          print(struc_break1)
        
      # break2
        returns_break2_1  =returns[time(returns) < struc_break1]
        returns_break2_2 =returns[time(returns) >=struc_break1]
        
        grid_struct_breaks2_1 = grid_struct_breaks[(grid_struct_breaks>=start_date+years(min_time_estim)) & (grid_struct_breaks <= (struc_break1-years(min_time_estim)))]
        grid_struct_breaks2_2 = grid_struct_breaks[(grid_struct_breaks<=end_date-years(min_time_estim)) & (grid_struct_breaks >= (struc_break1+years(min_time_estim)))]
        
        table_struc_break2_1 = find_structural_break(returns_break2_1,grid_struct_breaks2_1, start_parms, model_specification, number_restrictions, significance_level)
        table_struc_break2_2 = find_structural_break(returns_break2_2,grid_struct_breaks2_2, start_parms, model_specification, number_restrictions,significance_level)
        
        write.csv(table_struc_break2_1, file = paste0(outpathDescriptive,"table_struc_break2_1_",names(returns),".csv"))
        write.csv(table_struc_break2_2, file = paste0(outpathDescriptive,"table_struc_break2_2_",names(returns),".csv"))
        
      
        # manually select by looking at table
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
        
        table_struc_break3_1 = find_structural_break(returns_break3_1,grid_struct_breaks3_1, start_parms, model_specification, number_restrictions, significance_level)
        table_struc_break3_2 = find_structural_break(returns_break3_2,grid_struct_breaks3_2, start_parms, model_specification, number_restrictions,significance_level)
        table_struc_break3_3 = find_structural_break(returns_break3_3,grid_struct_breaks3_3, start_parms, model_specification, number_restrictions, significance_level)
        table_struc_break3_4 = find_structural_break(returns_break3_4,grid_struct_breaks3_4, start_parms, model_specification, number_restrictions,significance_level)
        
        write.csv(table_struc_break3_1, file = paste0(outpathDescriptive,"table_struc_break3_1_",names(returns),".csv"))
        write.csv(table_struc_break3_2, file = paste0(outpathDescriptive,"table_struc_break3_2_",names(returns),".csv"))
        write.csv(table_struc_break3_3, file = paste0(outpathDescriptive,"table_struc_break3_3_",names(returns),".csv"))
        write.csv(table_struc_break3_4, file = paste0(outpathDescriptive,"table_struc_break3_4_",names(returns),".csv"))
    }
  
  # save return for each timeperiod
    # only use struc_break1 and struc_break2_1 (struc_break2_2 is not sign for both. requires)
    ts_r_time1 = ts_r[(time(ts_r) < struc_break2_1)]
    ts_r_time2 = ts_r[(time(ts_r) >= struc_break2_1) & (time(ts_r)< struc_break1)]
    ts_r_time3 = ts_r[(time(ts_r) >= struc_break1)]
  
    all_garch_data_ts_r = list(ts_r_time1,ts_r_time2,ts_r_time3)
    names(all_garch_data_ts_r) = c("timeframe1", "timeframe2", "timeframe3")
    
    # option: use all structural breaks
      # ts_r_time3 = ts_r[(time(ts_r) >= struc_break1) & (time(ts_r)< struc_break2_2)]
      # ts_r_time4 = ts_r[(time(ts_r) >= struc_break2_2)]
      # all_garch_data_ts_r = list(ts_r_time1,ts_r_time2,ts_r_time3,ts_r_time4)
      # names(all_garch_data_ts_r) = c("timeframe1", "timeframe2", "timeframe3", "timeframe4")
    
# Autocorrelations of squares in timeframes ----
  # initalize lags and tables
    lags_list = c(1,5,10,20,40)
    lm_test_squares_struc_breaks = as.data.frame(matrix(data=NA,nrow=length(all_garch_data_ts_r)*2,ncol=(length(lags_list)+2)))
    colnames(lm_test_squares_struc_breaks) = c("Timeframe","Asset",paste0(lags_list, "lag(s)"))
    lm_test_squares_struc_breaks$Timeframe = c(rep(paste0(start_date,"/",(struc_break2_1-days(1))),2), 
                                               rep(paste0(struc_break2_1,"/",(struc_break1-days(1))),2),
                                               rep(paste0(struc_break1,"/",end_date),2))
    lm_test_squares_struc_breaks$Asset = c(rep(c("RUB/USD", "oil"),length(nrow(lm_test_squares_struc_breaks)/2)))
    lm_test_squares_struc_breaks_p_values = lm_test_squares_struc_breaks
  # test for each defined lag, asset and timeframe
  for (i in 1:length(lags_list)) {
    for(j in 1:length(all_garch_data_ts_r)) {
      #rub
      test_rub = ArchTest(all_garch_data_ts_r[[j]]$rub_errors, lags=lags_list[i], demean = TRUE)
      lm_test_squares_struc_breaks[(2*j-1),(i+2)] = test_rub$statistic
      lm_test_squares_struc_breaks_p_values[(2*j-1),(i+2)] = test_rub$p.value
      #oil
      test_oil = ArchTest(all_garch_data_ts_r[[j]]$oil_errors, lags=lags_list[i], demean = TRUE)
      lm_test_squares_struc_breaks[(2*j),(i+2)] = test_oil$statistic
      lm_test_squares_struc_breaks_p_values[(2*j),(i+2)] = test_oil$p.value
    }
  }
    
  # save test results
    write.csv(lm_test_squares_struc_breaks, file = paste0(outpathDescriptive,"lm_test_squares_struc_breaks",".csv"))
    write.csv(lm_test_squares_struc_breaks_p_values, file = paste0(outpathDescriptive,"lm_test_squares_struc_breaks_p_values",".csv"))
    xtable(lm_test_squares_struc_breaks, caption = "LM-test - teststatistic", digits=3)
    xtable(lm_test_squares_struc_breaks_p_values, caption = "LM-test - p-value", digits=3)
    
  # select timeframes for analysis: timeframe after struc_break1 
    garch_data_ts_r = all_garch_data_ts_r[[3]]
    saveRDS(garch_data_ts_r, file = paste0(outpathDescriptive,"garch_data_ts_r.rds"))
    
   print("Manually Select Timeframes with GARCH effects for analysis via LR-Test table (lm_test_squares_struc_breaks_p_values)")
   print("Selected Timeframe: 2008 until 2020, no strong structural breaks in this timeframe for both series (but indication)")
   

# Asymmetries for selected timeframe ----
  # test for asymmetries
      # plot and table for test of autocorrelation
      autocorr_asymmetries_oil = sampleAutocorrelation(garch_data_ts_r$oil_errors, "oil", 0.05, outpath= outpathDescriptive)
      autocorr_asymmetries_rub = sampleAutocorrelation(garch_data_ts_r$rub_errors, "RUBUSD", 0.05, outpath= outpathDescriptive)
      
      # save tables
       autocorr_asymmetries_oil_t = as.data.frame(t(as.matrix(autocorr_asymmetries_oil[,-1])))
        colnames(autocorr_asymmetries_oil_t) = paste0( colnames(autocorr_asymmetries_oil_t), "lag(s)")
        rownames(autocorr_asymmetries_oil_t) = c("corr epsilon+", "corr epsilon-", "p-val epsilon+", "p-val epsilon-")
  
        xtable(autocorr_asymmetries_oil_t, caption = "Asymmetries for oil - auto-correlation of positive and negative returns.", digits=3)
        write.csv(autocorr_asymmetries_oil_t, file = paste0(outpathDescriptive,"autocorr_asymmetries_oil_t",".csv"))
        
        autocorr_asymmetries_rub_t = as.data.frame(t(as.matrix(autocorr_asymmetries_rub[,-1])))
        colnames(autocorr_asymmetries_rub_t) =  paste0( colnames(autocorr_asymmetries_rub_t), "lag(s)")
        rownames(autocorr_asymmetries_rub_t) = c("corr epsilon+", "corr epsilon-", "p-val epsilon+", "p-val epsilon-")
        
        xtable(autocorr_asymmetries_rub_t, caption = "Asymmetries for RUB/USD - auto-correlation of positive and negative returns.", digits=3)
        write.csv(autocorr_asymmetries_rub_t, file = paste0(outpathDescriptive,"autocorr_asymmetries_rub_t",".csv"))

  # TODO! Asymmetries2: sign bias tests  ---- 

  # news impact curve? (see audrino code)