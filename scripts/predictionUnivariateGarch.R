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

# simple plots
ggplot(in_sample_pred_result$rub_tree,
       aes(x = date, y = variance_proxy)) +
  geom_line()
ggplot(in_sample_pred_result$rub_tree,
       aes(x = date, y = variance_predict)) +
  geom_line()
ggplot(in_sample_pred_result$rub_tree,
       aes(x = date, y = residuals_garch)) +
  geom_line()

 # geom_line(aes(y = variance_predict), color = "blue") 

ggplot(in_sample_pred_result$rub, aes(x = date, y = variance_proxy)) +
  geom_line(aes(y = variance_predict), color = "red") +
  geom_line(aes(y = residuals_garch), color = "blue")

ggplot(in_sample_pred_result$oil, aes(x = date, y = variance_proxy)) +
  geom_line(aes(y = variance_predict), color = "red") +
  geom_line(aes(y = residuals_garch), color = "blue")

# Correlations
cor(
  in_sample_pred_result$rub_tree$variance_predict,
  in_sample_pred_result$rub$variance_predict
)

cor(
  in_sample_pred_result$rub_tree$variance_predict,
  in_sample_pred_result$rub_tree$variance_proxy
)

cor(
  in_sample_pred_result$rub$variance_predict,
  in_sample_pred_result$rub$variance_proxy
)

cor(
  in_sample_pred_result$oil$variance_predict,
  in_sample_pred_result$oil$variance_proxy
)

# Residuals
  # squared residuals

  # distribution of residuals

