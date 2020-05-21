source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_v2.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_forecast.R")
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

####DCC Estimation Standard GARCH####
pred_results_1 <- full_sample[["predictions"]][["rub"]][["variance_predict"]]
pred_results_2 <- full_sample[["predictions"]][["oil"]][["variance_predict"]]
dccfit <- full_sample[["dccoutput"]][["fit"]]
returns <- full_sample[["dccoutput"]][["returns"]]
covs <- full_sample[["dccoutput"]][["fit"]]@mfit[["Q"]]

estimates <- dcc_forecast(dccfit, pred_results_1, pred_results_2, returns, covs)
full_sample$dccestimates <- estimates 



####TREE GARCH--------------------------------------------------------------------------
####Import Univariate GARCH Models####
rub_list_tree_1 <- full_sample[["models"]][["rub_tree"]][["rub_subsample1"]]
oil_list_tree_1 <- full_sample[["models"]][["oil"]][["oil_all"]]
rub_pred_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]]
oil_pred_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]]

####Call DCC Funcitons####
dccoutput_tree_1 <- dcc_function(rub_list_tree_1,oil_list_tree_1,rub_pred_tree_1,oil_pred_tree_1)
full_sample$dccoutput_tree_1 <- dccoutput_tree_1

####DCC Estimation Standard GARCH####
pred_results_1_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]][["variance_predict"]]
#pred_results_2_tree_1 <- full_sample[["predictions"]][["oil"]][["variance_predict"]] # i need them on the same date as rub
dccfit_tree_1 <- full_sample[["dccoutput"]][["fit"]]
returns_tree_1 <- full_sample[["dccoutput"]][["returns"]]
covs_tree_1 <- full_sample[["dccoutput"]][["fit"]]@mfit[["Q"]]

estimates_tree_1 <- dcc_forecast(dccfit_tree_1, pred_results_1_tree_1, pred_results_2_tree_1, returns_tree_1, covs_tree_1)
full_sample$dccestimates_tree_1 <- estimates_tree_1




####comparison####
returns <- full_sample[["dccoutput"]][["returns"]]
estimates_1 <- full_sample$dccestimates
estimates_2 <- full_sample$dccestimates
loss_function <- 3
output <- model_comparison(estimates_1, estimates_2, returns, loss_function = 3)

