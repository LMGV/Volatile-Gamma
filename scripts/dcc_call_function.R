source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_v2.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_forecast.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/model_comp.R")
library(zoo)
library(xts)
library(plyr)
####Import Univariate GARCH Models####
full_sample <- readRDS("C:/Users/johan/Documents/GitHub/Volatile-Gamma/output/univariateModels/model_and_prediction.rds")
full_sample[["predictions"]][["rub_tree_subsample2"]] <- full_sample[["predictions"]][["rub_tree_subsample2"]][-(1:3),]

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
####TREE Subsample 1####
####Import Univariate GARCH Models####
rub_list_tree_1 <- full_sample[["models"]][["rub_tree"]][["rub_subsample1"]]
oil_list_tree_1 <- full_sample[["models"]][["oil"]][["oil_all"]]
rub_pred_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]]
oil_pred_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]]

####Call DCC Funcitons####
dccoutput_tree_1 <- dcc_function(rub_list_tree_1,oil_list_tree_1,rub_pred_tree_1,oil_pred_tree_1)
full_sample$dccoutput_tree_1 <- dccoutput_tree_1

####DCC Estimation TREE-GARCH 1####
pred_results_1_tree_1 <- full_sample[["predictions"]][["rub_tree_subsample1"]][["variance_predict"]]
pred_results_2_tree_1 <- full_sample[["predictions"]][["oil"]]

ts<- xts(full_sample[["predictions"]][["rub_tree_subsample1"]], order.by = full_sample[["predictions"]][["rub_tree_subsample1"]]$date)
ts2 <- xts(full_sample[["predictions"]][["oil"]], order.by = full_sample[["predictions"]][["oil"]]$date)
 
pred_results_2_tree_1 <- merge(ts$date,ts2$variance_proxy,join='left')
pred_results_2_tree_1 <- as.numeric(pred_results_2_tree_1[,2])



dccfit_tree_1 <- full_sample[["dccoutput_tree_1"]][["fit"]]
returns_tree_1 <- full_sample[["dccoutput_tree_1"]][["returns"]]
covs_tree_1 <- full_sample[["dccoutput_tree_1"]][["fit"]]@mfit[["Q"]]

estimates_tree_1 <- dcc_forecast(dccfit_tree_1, pred_results_1_tree_1, pred_results_2_tree_1, returns_tree_1, covs_tree_1)
full_sample$dccestimates_tree_1 <- estimates_tree_1

####TREE Subsample 2####
####Import Univariate GARCH Models####
rub_list_tree_2 <- full_sample[["models"]][["rub_tree"]][["rub_subsample2"]]
oil_list_tree_2 <- full_sample[["models"]][["oil"]][["oil_all"]]
rub_pred_tree_2 <- full_sample[["predictions"]][["rub_tree_subsample2"]]
oil_pred_tree_2 <- full_sample[["predictions"]][["rub_tree_subsample2"]]

####Call DCC Funcitons####
dccoutput_tree_2 <- dcc_function(rub_list_tree_2,oil_list_tree_2,rub_pred_tree_2,oil_pred_tree_2)
full_sample$dccoutput_tree_2 <- dccoutput_tree_2

####DCC Estimation TREE-GARCH 2####
pred_results_1_tree_2 <- full_sample[["predictions"]][["rub_tree_subsample2"]][["variance_predict"]]
pred_results_2_tree_2 <- full_sample[["predictions"]][["oil"]]

ts <- xts(full_sample[["predictions"]][["rub_tree_subsample2"]], order.by = full_sample[["predictions"]][["rub_tree_subsample2"]]$date)
ts2 <- xts(full_sample[["predictions"]][["oil"]], order.by = full_sample[["predictions"]][["oil"]]$date)
pred_results_2_tree_2 <- merge(ts$date,ts2$variance_proxy,join='left')
pred_results_2_tree_2 <- as.numeric(pred_results_2_tree_2[,2])

dccfit_tree_2 <- full_sample[["dccoutput_tree_2"]][["fit"]]
returns_tree_2 <- full_sample[["dccoutput_tree_2"]][["returns"]]
covs_tree_2 <- full_sample[["dccoutput_tree_2"]][["fit"]]@mfit[["Q"]]

estimates_tree_2 <- dcc_forecast(dccfit_tree_2, pred_results_1_tree_2, pred_results_2_tree_2, returns_tree_2, covs_tree_2)
full_sample$dccestimates_tree_2 <- estimates_tree_2


####comparison####
est_1 <- as.data.frame(full_sample$dccestimates_tree_1)
rownames(est_1) <- full_sample[["predictions"]][["rub_tree_subsample1"]][["date"]]

est_2 <- as.data.frame(full_sample$dccestimates_tree_2)
rownames(est_2) <- full_sample[["predictions"]][["rub_tree_subsample2"]][["date"]]

est <- as.xts(rbind(est_1,est_2))
est_array <- array(data = , dim = c(nrow(est),2,2))
for (i in 1:nrow(est)) {
  est_array[i,1,1] <- est[i,1]
  est_array[i,2,1] <- est[i,2]
  est_array[i,1,2] <- est[i,3]
  est_array[i,2,2] <- est[i,4]
}




returns <- full_sample[["dccoutput"]][["returns"]]
estimates_1 <- full_sample$dccestimates
estimates_2 <- est_array
loss_function <- 3
output <- model_comparison(estimates_1, estimates_2, returns, loss_function = 3)

