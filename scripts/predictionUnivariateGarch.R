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
predict_data = data.frame(date = index(garch_data_ts_r), coredata(garch_data_ts_r)) %>%
  mutate(analysis_variable = NA)


model_in_sample_pred = list()
in_sample_pred_result = list()

model_list_choices = list(all_selected_model_tree, all_selected_model[1],all_selected_model[2] )
analysis_variable_choices = c("rub_errors", "oil_errors") # oil errors do not exist for tree so far

# 1) RUBUSE
# TREE MARCH model
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
  # !remove one missing obs in tree from other dataframes
  ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[2]]$date
  in_sample_pred_result[[2]] = in_sample_pred_result[[2]][-ind_same_date,]
  ind_same_date = in_sample_pred_result[[1]]$date %in% in_sample_pred_result[[3]]$date
  in_sample_pred_result[[3]] = in_sample_pred_result[[3]][-ind_same_date,]

# set names for series and save
names(in_sample_pred_result) = c("rub_tree", "rub", "oil")
names(model_in_sample_pred) = c("rub_tree", "rub", "oil")
saveRDS(in_sample_pred_result,
        file = paste0(outpathModels, "in_sample_pred_result.rds"))

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

# simple plots
ggplot(in_sample_pred_result$rub_tree,
       aes(x = date, y = variance_proxy)) +
  geom_line(aes(y = variance_predict), color = "red") +
  geom_line(aes(y = residuals_garch), color = "blue")

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

