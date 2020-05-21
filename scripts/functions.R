# basic functions

# in sample forecast given model, return data and dates
in_sample_forecast =  function(models, predict_data, start_date_predictions, end_date_predictions,max_lag_prediction) {
  predict_data_with_model = NULL
  model_list = list()
  for (subsample_iter in 1:length(models)){
    model_list[[subsample_iter]] = filter(predict_data, date %in% models[[subsample_iter]]$returns_with_date$date) %>%
      mutate(model = subsample_iter) #model number is number in original list
    
    predict_data_with_model = bind_rows(predict_data_with_model, model_list[[subsample_iter]]) # combine to one dataset again
  }
  
  # add volatility proxy - daily squared deviation
  predict_data_with_model= mutate(predict_data_with_model, variance_proxy = (analysis_variable-mean(rub,na.rm=T))^2,
                                  variance_predict = NA)
  
  # cut timeframe
  predict_data_with_model = filter(predict_data_with_model, date >=start_date_predictions & date <= end_date_predictions)
  
  # predict volatility for each observation 
  for (i in (2+max_lag_prediction):nrow(predict_data_with_model)){
    predict_data_with_model$variance_predict[i] = oneDayPredict(predict_data_with_model$analysis_variable[(i-max_lag_prediction-1):(i-1)], 
                                                                predict_data_with_model$variance_proxy[(i-max_lag_prediction-1):(i-1)],
                                                                models[[predict_data_with_model$model[i-1]]])
  }
  
  predict_data_with_model = predict_data_with_model[(max_lag_prediction+1):nrow(predict_data_with_model),] # remove missing obs due to lags
  predict_data_with_model$residuals_garch = predict_data_with_model$variance_predict- predict_data_with_model$variance_proxy
  results = select(predict_data_with_model, date, variance_proxy, variance_predict, residuals_garch, rub_errors, oil_errors,  rub_errors_lag1,  oil_errors_lag1)
  return(results)
  
}

# return volatility forecast for one day given model data. lag=1
oneDayPredict = function(past_returns, past_variance,model_data){
  # input: past max_lag obs of returns and variance, model data
  
  # get model specification and parms:
  model_coefs= model_data$garch_coefs
  model_specif= model_data$model_specification
  
  # input: for each day one return vector + list model coefs + list model spec  
  # output: scalar variance forecast
  
  mu_coef = model_coefs[1]
  constant_coef = model_coefs[2]
  ma_coef = model_coefs[(1:model_specif$number_ma)+2]
  ar_coef = model_coefs[((model_specif$number_ma+1):(model_specif$number_ma+model_specif$number_ar))+2]
  
  #threshhold for TGARCH
  if(model_specif$threshhold_included==T) {
    th_coef =  model_coefs$eta11 
    th_value = model_specif$th_value
  }
  
  # demeaned square return as variance proxy  
  epsilon  = past_returns -mu_coef$mu 
  epsilon_sq = epsilon^2
  
  # calc ma/ar/th parts of variance forecast
  ma_part = 0
  for (k in 1:length(ma_coef)) {
    ma_part = ma_part + ma_coef[k]*epsilon_sq[length(epsilon)+1-k]
  }
  ar_part = 0
  for (j in 1:length(ar_coef)) {
    ar_part = ar_part + ar_coef[j]*past_variance[length(epsilon)+1-j]  # calc AR part
  }
  
  # add treshhold parm if part of model
  threshhold_part = 0
  if(model_specif$threshhold_included ==T) {
    threshhold_part = th_coef*epsilon_sq[length(epsilon)]* as.numeric(epsilon[length(epsilon)]<=th_value)
  }
  # variance estimation
  var_estim = as.numeric(constant_coef$omega + ar_part+ ma_part+threshhold_part)
  
  return(var_estim)
  
}




## short helper functiions ----
my_aic = function(log_likelihood, number_parms) {
  aic = -2*log_likelihood+ 2*number_parms
  return(aic)
}

my_bic = function(log_likelihood, number_parms, T) {
  bic = -2*log_likelihood+ log(T)*number_parms
  return(bic)
}

value_at_risk_empirical = function(data, significance_level, time) {
  # data = returns in correct time scale
  # signif level= bound  for var, signifince
  # time = scaling factor for VaR with different timeframe. Requires normality
  if(missing(time)) { time=1}
  value_at_risk = -quantile(data,significance_level, na.rm=TRUE)*sqrt(time)
  return(value_at_risk)
}

expected_shortfall_empirical = function(data, significance_level, time) {
  # data = returns in correct time scale
  # signif level= bound  for var, signifince
  # time = scaling factor for VaR with different timeframe. Requires normality
  value_at_risk = -quantile(data,significance_level, na.rm=TRUE)
  expected_shortfall = -mean(subset(data, data < -value_at_risk), na.rm=TRUE)*sqrt(time)
  return(expected_shortfall)
}

supLikelihoodTest = function(loglik_unrestricted, loglik_restricted, number_restrictions, significanceLevel) { 
  # make sure input is loglikelihood, not negative loglikelihood
  p_value = 1-pchisq(2*(loglik_unrestricted-loglik_restricted),number_restrictions) 
  test_result = ifelse(p_value < significanceLevel, F,T)
  test_result_word = ifelse(p_value < significanceLevel, "h0_rejected","h0_not_rejected")
  result = c(test_result_word, test_result, p_value, 2*(loglik_unrestricted-loglik_restricted), loglik_unrestricted,loglik_restricted )
  names(result) = c("test_result_word","test_result","p_value","test_stat" , "logLikUR", "logLikR")
  return(result)
}

find_structural_break = function(returns,grid_struct_breaks, start_parms, model_specification, number_restrictions,significance_level) {
  # initialize test-result table
  struc_break_test_results = as.data.frame(matrix(nrow = length(grid_struct_breaks), ncol = 7,))
  colnames(struc_break_test_results) = c("break_date","test_result_word","test_result","p_value","test_stat" , "logLikUR", "logLikR")
  struc_break_test_results$break_date = grid_struct_breaks
  
  # estimate GARCH model for entire timeframe
  opt_parms= nlm(garchEstimation,start_parms,
                 returns = returns,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                 threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                 distribution=model_specification$distribution,
                 print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
  
  # estimate 2 GARCH models for before and after structural breaks in grid 
  for (break_iter in 1:length(grid_struct_breaks)) {
    sample_before = returns[index(returns) < grid_struct_breaks[break_iter]]
    sample_after   = returns[index(returns) >= grid_struct_breaks[break_iter]]
    
    # estimate GARCH model before and after potential breakpoint
    opt_parms_before= nlm(garchEstimation,start_parms,
                          returns = sample_before,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                          threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                          distribution=model_specification$distribution,
                          print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
    
    opt_parms_after= nlm(garchEstimation,start_parms,
                         returns = sample_after,  ma = model_specification$number_ma, ar = model_specification$number_ar, 
                         threshhold = model_specification$threshhold_included, th_value = model_specification$th_value, data_threshhold = data_threshhold,
                         distribution=model_specification$distribution,
                         print.level=0,steptol = 1e-6, iterlim=1000, check.analyticals=T)
    
    # LR-test
    struc_break_test_results[break_iter,2:7] = supLikelihoodTest(-opt_parms_before$minimum - opt_parms_after$minimum,-opt_parms$minimum, number_restrictions,significance_level) 
  }
  return(struc_break_test_results)
}   

## sample autocorrelation plots / table
sampleAutocorrelation = function(returns, asset_name, significance_level, outpath) {
  epsilon_pos = pmax((returns-mean(returns)),0)
  epsilon_neg = pmin((returns-mean(returns)),0)
  
  # sample autocorrelation
  list_lags = c(1,2,5,10,20,40) 
  num_lags = max(list_lags)
  
  lag = seq(1, num_lags)
  corr_pos = rep(0, num_lags)
  corr_neg = rep(0, num_lags)
  corr = as.data.frame(cbind(lag,corr_pos,corr_neg))
  for (h in 1:num_lags)
  {
    corr[h,2] =  cor(epsilon_pos[(1+h):(length(returns))],(returns-mean(returns))[(1):(length(returns)-h)]) #Leverage effect
    corr[h,3] = cor(epsilon_neg[(1+h):(length(returns))],(returns-mean(returns))[(1):(length(returns)-h)]) #"reverse" leverage effect
  }
  
  # significance autocorrelation
  # approximate variance with 1/T, 2sided test
  corr$corr_pos_p_value  = 1-pnorm(abs(corr$corr_pos *sqrt(length(corr$corr_pos)) * sqrt(length(corr$corr_pos)))) # mu * sqrt(N) * 1/sigma(corr_pos)
  corr$corr_neg_p_value  = 1-pnorm(abs(corr$corr_neg *sqrt(length(corr$corr_neg)) * sqrt(length(corr$corr_neg)))) # mu * sqrt(N) * 1/sigma(corr_neg)
  
  title= paste0("Asymmetries_", asset_name)
  x = corr$lag
  xlab  ="lags"
  ylab = "correlation"
  names_y = colnames(corr)[1:2]
  y1 = corr$corr_pos
  y2 = corr$corr_neg
  legend = T
  y_percent = F
  y_discrete = F
  print(line_plot_multiple(title, outpath,x,xlab, ylab, names_y, y_percent, y_discrete, legend, y1, y2))
  
  title= paste0("Asymmetries_P_values ", asset_name)
  x = corr$lag
  xlab  ="lags"
  ylab = "p_value"
  names_y = c(colnames(corr)[3:4], "sign_bound")
  y1 = corr$corr_pos_p_value
  y2 = corr$corr_neg_p_value
  y3 = significance_level/2
  legend = T
  y_percent = F
  y_discrete = F
  line_plot_multiple(title, outpath, x,xlab, ylab, names_y, y_percent, y_discrete, legend, y1, y2,y3)
  
  # prepare table
  # only include certain lags in table
  corr_table = corr[list_lags,]
  return(corr_table)
}

## plotting function ----
line_plot_multiple = function(title, outpath, x,xlab, ylab, names_y, y_percent, y_discrete, legend, y1, y2,y3,y4,y5,y6,y7,y8,y9,y10) {
  # plots up to 10 lines with same scale
  # insert vectors for x, y1, y2... y2-10 are optional
  # give vector for legend entry names_y, ...
  # if names_y na or missing colors are still given
  # legend: boolean if legend should be created, Default is TRUE
  if(missing(xlab)) {xlab=NULL}
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  if(missing(names_y)==T | all(is.na(names_y))==T) {names_y=c(1:10)}
  if(missing(legend)){legend = T}
  if(missing(y_percent)){y_percent = FALSE}

  plot = ggplot(data=NULL, aes(x=x))+
    geom_line(aes(y=y1, col = names_y[1]))

  # only add layers if value provided
  if (missing(y2)==F) {
    plot = plot + geom_line(aes(y =y2 , col =names_y[2]))
  }
  if (missing(y3)==F) {
    plot = plot + geom_line(aes(y =y3 , col =names_y[3]))
  }
  if (missing(y4)==F) {
    plot = plot + geom_line(aes(y =y4 , col =names_y[4]))
  }
  if (missing(y5)==F) {
    plot = plot + geom_line(aes(y =y5 , col =names_y[5]))
  }
  if (missing(y6)==F) {
    plot = plot + geom_line(aes(y =y6 , col =names_y[6]))
  }
  if (missing(y7)==F) {
    plot = plot + geom_line(aes(y =y7 , col =names_y[7]))
  }
  if (missing(y8)==F) {
    plot = plot + geom_line(aes(y =y8 , col =names_y[8]))
  }
  if (missing(y9)==F) {
    plot = plot + geom_line(aes(y =y9 , col =names_y[9]))
  }
  if (missing(y10)==F) {
    plot = plot + geom_line(aes(y =y10 , col =names_y[10]))
  }


  plot = plot + theme_bw()+
    labs(x = xlab)+
    labs(y = ylab)

  # add legend if Legend is true are given
  if(legend==T) {
    plot = plot + theme(legend.title = element_blank(), legend.position = "bottom", legend.box.background = element_rect(colour = "black"))
  } else {
    plot = plot + theme(legend.position = "none")
  }

  # add units to yaxis (e.g. percent)
  if(y_percent==T) {
    plot = plot + scale_y_continuous(labels = scales::percent_format(accuracy = 2))
  }

  if(y_discrete==T) {
    plot = plot + scale_y_discrete()
  }

  # other plot options
  plot = plot +  ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0(outpath,title,".png"), width=6, height=4, dpi=600)


  print(plot)
  return(plot)
}

line_point_plot_multiple = function(title, x,xlab, ylab, names_y, y_percent,  y_discrete, legend, y1, y2,y3,y4,y5,y6,y7,y8,y9,y10) {
  # plots up to 10 lines with same scale 
  # insert vectors for x, y1, y2... y2-10 are optional
  # give vector for legend entry names_y, ...
  # if names_y na or missing colors are still given
  # legend: boolean if legend should be created, Default is TRUE
  # y6-y10 are geom point
  if(missing(xlab)) {xlab=NULL}
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  if(missing(names_y)==T | all(is.na(names_y))==T) {names_y=c(1:10)}
  if(missing(legend)){legend = T}
  if(missing(y_percent)){y_percent = FALSE}
  
  #first geom point
  plot = ggplot(data=NULL, aes(x=x))+
    geom_point(aes(y=y1, col = names_y[1])) 
  
  # geom lines
  if (missing(y2)==F) {
    plot = plot + geom_line(aes(y =y2 , col =names_y[2]))
  }  
  if (missing(y3)==F) {
    plot = plot + geom_line(aes(y =y3 , col =names_y[3]))
  }  
  if (missing(y4)==F) {
    plot = plot + geom_line(aes(y =y4 , col =names_y[4]))
  }  
  if (missing(y5)==F) {
    plot = plot + geom_line(aes(y =y5 , col =names_y[5]))
  }  
  
  # geom points
  if (missing(y6)==F) {
    plot = plot + geom_point(aes(y =y6 , col =names_y[6]))
  }  
  if (missing(y7)==F) {
    plot = plot + geom_point(aes(y =y7 , col =names_y[7]))
  }  
  if (missing(y8)==F) {
    plot = plot + geom_point(aes(y =y8 , col =names_y[8]))
  }  
  if (missing(y9)==F) {
    plot = plot + geom_point(aes(y =y9 , col =names_y[9]))
  }  
  if (missing(y10)==F) {
    plot = plot + geom_point(aes(y =y10 , col =names_y[10]))
  }  
  
  
  plot = plot + theme_bw()+
    labs(x = xlab)+
    labs(y = ylab)
  
  # add legend if Legend is true are given
  if(legend==T) {
    plot = plot + theme(legend.title = element_blank(), legend.position = "bottom", legend.box.background = element_rect(colour = "black"))
  } else {
    plot = plot + theme(legend.position = "none")
  }
  
  # add units to yaxis (e.g. percent)
  if(y_percent==T) {
    plot = plot + scale_y_continuous(labels = scales::percent_format(accuracy = 2))
  }
  
  # y discrete=
  if(y_discrete==T) {
    plot = plot + scale_y_discrete()
  }
  
  # other plot options
  plot = plot +  ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+ 
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  
  print(plot)
  return(plot)
}


# times series plot function
time_series_plot_basic = function(date,title, y, ylab, min_date, max_date) {
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  
  plot = ggplot(data=NULL, aes(x=as.Date(date,format = "%m/%d/%Y"), y = y))+
    geom_line(colour = "cornflowerblue", size = 1)+
    scale_x_date(limits = c(min=min_date, max=max_date)) + 
    theme_bw()+
    labs(x = "Date")+
    labs(y = ylab)+
    theme(legend.position = "bottom", legend.box.background = element_rect(colour = "black"))+
    ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  plot
  return(plot)
}



time_series_plot_multiple = function(date,title, y, ylab, seperating_var,min_date, max_date) {
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  
  plot = ggplot(data=NULL, aes(x=as.Date(date,format = "%m/%d/%Y"), y = y))+
    geom_line(aes(color = seperating_var), size = 1)+
    scale_x_date(limits = c(min=min_date, max=max_date)) + 
    theme_bw()+
    labs(x = "Date")+
    labs(y = ylab)+
    theme(legend.position = "bottom", legend.box.background = element_rect(colour = "black"))+
    ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  
  print(plot)
  return(plot)
}
