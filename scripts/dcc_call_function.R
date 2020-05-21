source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_v2.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_forecasts.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/model_comp.R")


####Import Univariate GARCH Models####
full_sample <- readRDS("C:/Users/johan/Documents/GitHub/Volatile-Gamma/output/univariateModels/model_and_prediction.rds")

rub_list <- full_sample[["models"]][["rub"]][["rub_all"]]
oil_list <- full_sample[["models"]][["oil"]][["oil_all"]]
rub_pred <- full_sample[["predictions"]][["rub"]]
oil_pred <- full_sample[["predictions"]][["oil"]]

####Call DCC Funcitons####
dccoutput <- dcc_function(rub_list,oil_list,rub_pred,oil_pred)
full_sample$dccoutput <- dccoutput

####DCC Estimation
pred_results_1 <- full_sample[["predictions"]][["rub"]][["variance_predict"]]
pred_results_2 <- full_sample[["predictions"]][["oil"]][["variance_predict"]]
dccfit <- full_sample[["dccoutput"]][["fit"]]
returns <- full_sample[["dccoutput"]][["returns"]]
covs <- full_sample[["dccoutput"]][["fit"]]@mfit[["Q"]]

estimates <- dcc_forecast(dccfit, pred_results_1, pred_results_2, returns, covs)
full_sample$dccestimates <- estimates 

####comparison####
returns <- full_sample[["dccoutput"]][["returns"]]
estimates_1 <- full_sample$dccestimates
estimates_2 <- full_sample$dccestimates
loss_function <- 3
output <- model_comparison(estimates_1, estimates_2, returns, loss_function = 3)

