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
library(car)
library(sandwich)
library(lmtest)

filter <- dplyr::filter
select <- dplyr::select
setwd("~/GitHub/Volatile-Gamma") # setwd

source("scripts/functions.R") # functions
source("scripts/garchFunction.R") # functions
outpathDescriptive = "output/univariateDescriptives/"
outpathModels =  "output/univariateModels/"

# Import Data ----
garch_data_ts_r  = readRDS("output/univariateDescriptives/garch_data_ts_r.rds") # selected garch data after struc break analysis
garch_data_ts_r_errors = garch_data_ts_r[, c("rub_errors", "oil_errors")]

# load model data:
  # static models
  all_selected_model  = readRDS("output/univariateModels/univariate_garchs_full_sample.rds")
  all_selected_model_tree  = readRDS("output/univariateModels/univariate_garchs_tree_rub.rds")
  all_selected_model_custom = readRDS("output/univariateModels/univariate_garchs_custom_rub.rds")

# DEFINE split_variables
# only lags etc that are included
# garch_data_ts_r should contain all vars
# garch_data_ts_r_errors only returns itself

# Prediction ----
# DEFINE split_variables
# only lags etc that are included
# garch_data_ts_r should contain all vars
# garch_data_ts_r_errors only returns itself



# in sample 1d-head forecasting for RUB
max_lag_prediction = 3
start_date_predictions = "2008-01-10"
end_date_predictions = "2019-12-31"
analysis_variable = c("rub_errors", "oil_errors")


# tree GARCH
# assign models to observations
garch_data_ts_r$rub_errors_lag1 = lag(garch_data_ts_r$rub_errors,1)
garch_data_ts_r$oil_errors_lag1 = lag(garch_data_ts_r$oil_errors,1)

predict_data = data.frame(date = index(garch_data_ts_r), coredata(garch_data_ts_r)) %>%
  mutate(analysis_variable = NA)


model_in_sample_pred = list()
in_sample_pred_result = list()

model_in_sample_pred_model_compare = list()
in_sample_pred_result_model_compare = list()

model_list_choices = list(all_selected_model_tree,all_selected_model[1],all_selected_model[2])
model_compare_list_choices = list(all_selected_model_custom[1],all_selected_model_custom[2],all_selected_model_custom[3])
analysis_variable_choices = c("rub_errors", "oil_errors") # oil errors do not exist for tree so far

# 1) RUBUSE
# TREE GARCH model
# define inputs for function
models = model_list_choices[[1]]
analysis_variable = analysis_variable_choices[1]
predict_data[, colnames(predict_data) == "analysis_variable"] = predict_data[, colnames(predict_data) == analysis_variable]

# predict variance
model_in_sample_pred[[1]] = models
in_sample_pred_result[[1]] = in_sample_forecast(
  models = models,
  predict_data,
  start_date_predictions,
  end_date_predictions,
  max_lag_prediction
)


# FULL Sample GARCH model
# define inputs for function
models = model_list_choices[[2]]
analysis_variable = analysis_variable_choices[1]
predict_data[, colnames(predict_data) == "analysis_variable"] = predict_data[, colnames(predict_data) == analysis_variable]

# predict variance
model_in_sample_pred[[2]] = models
in_sample_pred_result[[2]] = in_sample_forecast(
  models = models ,
  predict_data,
  start_date_predictions,
  end_date_predictions,
  max_lag_prediction
)

# Full Sample Different models rub
  for (i in 1:length(model_compare_list_choices)) {
  # model1
  models = model_compare_list_choices[[i]]
  analysis_variable = analysis_variable_choices[1] #same analysis variable
  predict_data[, colnames(predict_data) == "analysis_variable"] = predict_data[, colnames(predict_data) == analysis_variable]
  
  # predict variance
  model_in_sample_pred_model_compare[[i]] = models
  in_sample_pred_result_model_compare[[i]] = in_sample_forecast(
    models = models ,
    predict_data,
    start_date_predictions,
    end_date_predictions,
    max_lag_prediction
  )
}


# 2) OIL
# FULL Sample GARCH model
# define inputs for function
models = model_list_choices[[3]]
analysis_variable = analysis_variable_choices[2]
predict_data[, colnames(predict_data) == "analysis_variable"] = predict_data[, colnames(predict_data) == analysis_variable]

# predict variance
model_in_sample_pred[[3]] = models
in_sample_pred_result[[3]] = in_sample_forecast(
  models = models ,
  predict_data,
  start_date_predictions,
  end_date_predictions,
  max_lag_prediction
)

# match dates in in_sample_pred_result
  # remove missing forecasts
  in_sample_pred_result[[1]] = in_sample_pred_result[[1]][is.na(in_sample_pred_result[[1]]$variance_predict)==F,]
  in_sample_pred_result[[2]] = in_sample_pred_result[[2]][is.na(in_sample_pred_result[[2]]$variance_predict)==F,]
  in_sample_pred_result[[3]] = in_sample_pred_result[[3]][is.na(in_sample_pred_result[[3]]$variance_predict)==F,]

  # !remove one missing obs in tree from other dataframes
  ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[2]]$date
  in_sample_pred_result[[2]] = in_sample_pred_result[[2]][-ind_same_date,]
  ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[3]]$date
  in_sample_pred_result[[3]] = in_sample_pred_result[[3]][-ind_same_date,]
  
  ## EDIT FOR DCC input
  # add seperate predictions for subsapmle 1 and 2 of tree and remove missings
    rub_tree_subsample1 = left_join(all_selected_model_tree$rub_subsample1$returns_with_date, in_sample_pred_result[[1]], by ="date")
    rub_tree_subsample1 = rub_tree_subsample1[is.na(rub_tree_subsample1$variance_predict)==F,]
    rub_tree_subsample2 = left_join(all_selected_model_tree$rub_subsample2$returns_with_date, in_sample_pred_result[[1]], by ="date")
    rub_tree_subsample2 = rub_tree_subsample2[is.na(rub_tree_subsample2$variance_predict)==F,]
  
    in_sample_pred_result[[4]] = rub_tree_subsample1
    in_sample_pred_result[[5]] = rub_tree_subsample2
    
# match dates in in_sample_pred_result
# remove missing forecasts
in_sample_pred_result[[1]] = in_sample_pred_result[[1]][is.na(in_sample_pred_result[[1]]$variance_predict)==F,]
in_sample_pred_result[[2]] = in_sample_pred_result[[2]][is.na(in_sample_pred_result[[2]]$variance_predict)==F,]
in_sample_pred_result[[3]] = in_sample_pred_result[[3]][is.na(in_sample_pred_result[[3]]$variance_predict)==F,]

# !remove one missing obs in tree from other dataframes
ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[2]]$date
in_sample_pred_result[[2]] = in_sample_pred_result[[2]][-ind_same_date,]
ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[3]]$date
in_sample_pred_result[[3]] = in_sample_pred_result[[3]][-ind_same_date,]


## univ model compare: add tree forecasts
  in_sample_pred_result_model_compare[[4]] = in_sample_pred_result[[1]]

  # univ model compare: remove missing forecasts
  in_sample_pred_result_model_compare[[1]] = in_sample_pred_result_model_compare[[1]][is.na(in_sample_pred_result_model_compare[[1]]$variance_predict)==F,]
  in_sample_pred_result_model_compare[[2]] = in_sample_pred_result_model_compare[[2]][is.na(in_sample_pred_result_model_compare[[2]]$variance_predict)==F,]
  in_sample_pred_result_model_compare[[3]] = in_sample_pred_result_model_compare[[3]][is.na(in_sample_pred_result_model_compare[[3]]$variance_predict)==F,]

  # match dates of tree
  ind_same_date = in_sample_pred_result_model_compare[[4]]$date %in% in_sample_pred_result_model_compare[[1]]$date
  in_sample_pred_result_model_compare[[1]] = in_sample_pred_result_model_compare[[1]][-ind_same_date,]
  ind_same_date = in_sample_pred_result_model_compare[[4]]$date %in% in_sample_pred_result_model_compare[[2]]$date
  in_sample_pred_result_model_compare[[2]] = in_sample_pred_result_model_compare[[2]][-ind_same_date,]
  ind_same_date = in_sample_pred_result_model_compare[[4]]$date %in% in_sample_pred_result_model_compare[[3]]$date
  in_sample_pred_result_model_compare[[3]] = in_sample_pred_result_model_compare[[3]][-ind_same_date,]
  ind_same_date = in_sample_pred_result_model_compare[[4]]$date %in% in_sample_pred_result_model_compare[[4]]$date
  in_sample_pred_result_model_compare[[4]] = in_sample_pred_result_model_compare[[4]][-ind_same_date,]
  
  
# set names for series and save
#for dcc
names(in_sample_pred_result) = c("rub_tree", "rub", "oil","rub_tree_subsample1","rub_tree_subsample2")
names(model_in_sample_pred) = c("rub_tree", "rub", "oil")
saveRDS(in_sample_pred_result,
        file = paste0(outpathModels, "in_sample_pred_result.rds"))

# for univariate
names(in_sample_pred_result_model_compare) = c("rub_normal_11","rub_t_11" , "gjr_rub_t_11","rub_tree")
names(model_in_sample_pred_model_compare) = c("rub_normal_11","rub_t_11" , "gjr_rub_t_11")
saveRDS(in_sample_pred_result,
        file = paste0(outpathModels, "in_sample_pred_result_model_compare.rds"))


# list with models and predictions
  overall_model_list = list(model_in_sample_pred, in_sample_pred_result)
  names(overall_model_list) = c("models","predictions")
  saveRDS(overall_model_list,
        file = paste0(outpathModels, "model_and_prediction.rds"))

summary(
  in_sample_pred_result$rub_tree$variance_predict - in_sample_pred_result$rub$variance_predict
)

# brief check of results
for (i in 1:length(in_sample_pred_result)) {
  print(names(in_sample_pred_result)[i])
  print(summary(in_sample_pred_result[[i]]$variance_predict))
  print(summary(in_sample_pred_result[[i]]$residuals_garch))
}


# out of sample forecasting - not done ---
# estimate model each quarter
  estimation_timeframe = 1 # length estimation time in years
  grid_estimation_time = seq(min(time(garch_data_ts_r_errors)) + years(estimation_timeframe), max(time(garch_data_ts_r_errors)), by = "quarter")
  
  # estimate model each quarter
  for (quarter in grid_estimation_time) {
    
  }


# Prediction Analysis ----
  summary(in_sample_pred_result$rub_tree$variance_proxy)
  summary(in_sample_pred_result$rub_tree$variance_proxy)
  summary(in_sample_pred_result$rub_tree$residuals_garch)
  summary(in_sample_pred_result$rub$residuals_garch)
  summary(in_sample_pred_result$oil$residuals_garch)

# Plot Vol proxys and model ----
  
  # add MAE proxy
  for (i in 1:length(in_sample_pred_result_model_compare))  {
    in_sample_pred_result_model_compare[[i]]$mae_proxy = abs(in_sample_pred_result_model_compare[[i]]$rub_errors)
    in_sample_pred_result_model_compare[[i]]$vol_proxy_standardized = in_sample_pred_result_model_compare[[i]]$variance_proxy /in_sample_pred_result_model_compare[[i]]$variance_predict
    in_sample_pred_result_model_compare[[i]]$mae_proxy_standardized = in_sample_pred_result_model_compare[[i]]$mae_proxy /in_sample_pred_result_model_compare[[i]]$variance_predict
    }
  # plot all returns
  
  title = "RUBUSD Daily Returns Volatility"
  xlab = "Time"
  ylab= ""
  y1 = sqrt(in_sample_pred_result_model_compare[[1]]$variance_proxy)
  x = in_sample_pred_result_model_compare[[1]]$date
  names_y = c("RUB/USD")
  line_plot_multiple(title, outpath = outpathModels, x,xlab, ylab, names_y, y_percent=T, y_discrete=F, legend=F, y1)
  
  title = "RUBUSD Abs Daily Returns"
  xlab = "Time"
  ylab= ""
  y1 = in_sample_pred_result_model_compare[[1]]$mae_proxy
  x = in_sample_pred_result_model_compare[[1]]$date
  names_y = c("RUB/USD")
  line_plot_multiple(title, outpath = outpathModels, x,xlab, ylab, names_y, y_percent=T, y_discrete=F, legend=F, y1)
  
  title_plot = paste("RUBUSD Volatility -", c("Garch(1,1) normal", "Garch(1,1) t", "GJR-Garch(1,1) t", "Tree-Garch(1,1) t"))
  for (i in 1:length(in_sample_pred_result_model_compare))  {
    
    title = title_plot[i]
    xlab = "Time"
    ylab= ""
    y1 = sqrt(in_sample_pred_result_model_compare[[i]]$variance_predict)
    x = in_sample_pred_result_model_compare[[i]]$date
    names_y = c("RUB/USD")
    line_plot_multiple(title, outpath = outpathModels, x,xlab, ylab, names_y, y_percent=T, y_discrete=F, legend=F, y1)
    
  }


# Minzer-Zarnowitz-regressions for RUB ----
data_predictions = in_sample_pred_result_model_compare[-5]

test_results = as.data.frame(matrix(nrow= 4, ncol = 4,))
colnames(test_results) = c("model","alpha","beta","R^2 MZ")
test_results$model = c("Garch(1,1) normal", "Garch(1,1) t", "GJR-Garch(1,1) t", "Tree-Garch(1,1) t")

test_results_pval = as.data.frame(matrix(nrow= 4, ncol = 4,))
colnames(test_results_pval) = c("model","alpha","beta", "joint_test")
test_results_pval$model = c("Garch(1,1) normal", "Garch(1,1) t", "GJR-Garch(1,1) t", "Tree-Garch(1,1) t")

  for (i in 1:length(data_predictions)) {
data_predictions_selected = data_predictions[[i]]

  # classic regression
  classic_mz = lm(variance_proxy ~ variance_predict, data = data_predictions_selected)
  test_results[i,2:3] = coeftest(classic_mz, vcov = vcovHC(classic_mz, type="HC1"))[1:2,1]
  test_results[i,4] = summary(classic_mz)$r.squared
  test_results_pval[i,2] = car::linearHypothesis(classic_mz, c("(Intercept) = 0"))$`Pr(>F)`[2]
  test_results_pval[i,3] = car::linearHypothesis(classic_mz, c("variance_predict = 1"))$`Pr(>F)`[2]
  test_results_pval[i,4] = car::linearHypothesis(classic_mz, c("(Intercept) = 0" , "variance_predict = 1"))$`Pr(>F)`[2]
  
  }


xtable(test_results, caption = "MZ - Regression with squared return proxy", digits = 3)
xtable(test_results_pval, caption = "PVALUES for MZ - Regression with squared return proxy", digits = 3)
  

# Likelihood tables

lik_tree = all_selected_model_tree$rub_subsample1$model_evaluation$log_lik + all_selected_model_tree$rub_subsample2$model_evaluation$log_lik
aic_tree = all_selected_model_tree$rub_subsample1$model_evaluation$aic_model + all_selected_model_tree$rub_subsample2$model_evaluation$aic_model

likelihood_models = as.data.frame(matrix(nrow=2, ncol=5))
colnames(likelihood_models) = c("","Garch(1,1) normal", "Garch(1,1) t", "GJR-Garch(1,1) t", "Tree-Garch(1,1) t")
likelihood_models[,1] = c("logLik","AIC")


likelihood_models$`Garch(1,1) normal` = c(all_selected_model_custom$rub_normal_11$model_evaluation$log_lik, all_selected_model_custom$rub_normal_11$model_evaluation$aic_model)
likelihood_models$`Garch(1,1) t` = c(all_selected_model_custom$rub_t_11$model_evaluation$log_lik, all_selected_model_custom$rub_t_11$model_evaluation$aic_model)
  likelihood_models$`GJR-Garch(1,1) t` =  c(all_selected_model_custom$rub_t_gjr_11$model_evaluation$log_lik, all_selected_model_custom$rub_t_gjr_11$model_evaluation$aic_model)
  likelihood_models$`Tree-Garch(1,1) t`= c(lik_tree,aic_tree)
  
xtable(likelihood_models, caption = "Likelihood comparision univariate GARCH-models for RUBUSD", ditigs=3)

# Parameters Garchs

coef_table = as.data.frame(matrix(nrow=5, ncol=4))
colnames(coef_table) = c("model","alpha","beta","df_t")
coef_table$model = c("Garch(1,1) normal", "Garch(1,1) t", "GJR-Garch(1,1) t", "Subsample 1 -Tree-Garch", "Subsample 2 -Tree-Garch")

  coef_table[1,2:4] = c(all_selected_model_custom$rub_normal_11$garch_coefs[c("alpha1","beta1")],"/")
  coef_table[2,2:4] = all_selected_model_custom$rub_t_11$garch_coefs[c("alpha1","beta1","shape")]
  coef_table[3,2:4] = all_selected_model_custom$rub_t_gjr_11$garch_coefs[c("alpha1","beta1","shape")]
  coef_table[4,2:4] = all_selected_model_tree$rub_subsample1$garch_coefs[c("alpha1","beta1","shape")]
  coef_table[5,2:4] = all_selected_model_tree$rub_subsample2$garch_coefs[c("alpha1","beta1","shape")]

  xtable(coef_table, caption = "Coefficient comparision univariate GARCH-models for RUBUSD",digits =3)
  
  