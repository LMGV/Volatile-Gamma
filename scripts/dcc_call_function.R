source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_v2.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/dcc_forecast.R")
source("C:/Users/johan/Documents/GitHub/Volatile-Gamma/scripts/model_comp.R")
library(zoo)
library(xts)
library(plyr)
library(dplyr)
library(tidyr)
library(xtable)

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
ts_tree <- index(est)
est_array <- array(data = 0, dim = c(nrow(est),2,2))
for (i in 1:nrow(est)) {
  est_array[i,1,1] <- est[i,1]
  est_array[i,2,1] <- est[i,2]
  est_array[i,1,2] <- est[i,3]
  est_array[i,2,2] <- est[i,4]
}

estimates_1 <- full_sample$dccestimates
estimates_2 <- est_array

estimates_1 <- as.data.frame(estimates_1)
estimates_2 <- as.data.frame(estimates_2)
returns <- full_sample[["dccoutput"]][["returns"]]
ts <- full_sample[["predictions"]][["rub"]][["date"]]

estimates_1 <- xts(estimates_1, order.by = ts)
estimates_11 <- cbind(index(estimates_1), as.data.frame(estimates_1))
colnames(estimates_11)[1] <- "Date" 
estimates_111 <- estimates_11 %>% separate(Date, c("Date", "tzone"), " ")

estimates_2 <- xts(estimates_2, order.by = index(est))
estimates_22 <- cbind(index(estimates_2), as.data.frame(estimates_2))
colnames(estimates_22)[1] <- "Date" 
estimates_222 <- estimates_22 %>% separate(Date, c("Date", "tzone"), " ")



returns <- xts(returns, order.by = ts)
returns_1 <- cbind(index(returns), as.data.frame(returns))
colnames(returns_1)[1] <- "Date" 
returns_11 <- returns_1 %>% separate(Date, c("Date", "tzone"), " ")


merged_estimates <- inner_join(estimates_111, estimates_222, by="Date")
merged_estimates_returns <- inner_join(merged_estimates, returns_11, by="Date")

estimates_1 <- array(data = 0, dim = c(nrow(merged_estimates),2,2))
estimates_2 <- array(data = 0, dim = c(nrow(merged_estimates),2,2))
for (i in 1:nrow(est)) {
  estimates_1[i,1,1] <- merged_estimates_returns[i,3]
  estimates_1[i,2,1] <- merged_estimates_returns[i,4]
  estimates_1[i,1,2] <- merged_estimates_returns[i,5]
  estimates_1[i,2,2] <- merged_estimates_returns[i,6]
  estimates_2[i,1,1] <- merged_estimates_returns[i,8]
  estimates_2[i,2,1] <- merged_estimates_returns[i,9]
  estimates_2[i,1,2] <- merged_estimates_returns[i,10]
  estimates_2[i,2,2] <- merged_estimates_returns[i,11]
}

returns <- merged_estimates_returns[,13:14]

output_MSE <- model_comparison(estimates_1, estimates_2, returns, loss_function = 3)
output_MSE

output_MEA <- model_comparison(estimates_1, estimates_2, returns, loss_function = 2)

####Table output####
#DMW
model_comp <- matrix(NA, nrow = 2, ncol = 4)
model_comp[1:2,1:2] <- output_MSE$Loss_function
model_comp[1,3:4] <- output_MSE$DMW
model_comp[2,3:4] <- output_MEA$DMW
colnames(model_comp) <- c("GARCH", "Tree-GARCH", "DMW", "p")
rownames(model_comp) <- c("MSE", "MEA")
xtable(model_comp,digits = 4)

#DCC
dcc_pars <- matrix(NA, nrow = 3, ncol = 2)
dcc_pars[1,1] <- full_sample[["dccoutput"]][["fit"]]@mfit[["coef"]][["[Joint]dcca1"]]
dcc_pars[2,1] <- full_sample[["dccoutput_tree_1"]][["fit"]]@mfit[["coef"]][["[Joint]dcca1"]]
dcc_pars[3,1] <- full_sample[["dccoutput_tree_2"]][["fit"]]@mfit[["coef"]][["[Joint]dcca1"]]

dcc_pars[1,2] <- full_sample[["dccoutput"]][["fit"]]@mfit[["coef"]][["[Joint]dccb1"]]
dcc_pars[2,2] <- full_sample[["dccoutput_tree_1"]][["fit"]]@mfit[["coef"]][["[Joint]dccb1"]]
dcc_pars[3,2] <- full_sample[["dccoutput_tree_2"]][["fit"]]@mfit[["coef"]][["[Joint]dccb1"]]

colnames(dcc_pars) <- c("alpha1","beta1")
rownames(dcc_pars) <- c("dcc-GARCH", "dcc-Tree-GARCH sample 1","dcc-Tree-GARCH sample 2")
xtable(dcc_pars,digits = 5)

#Info criterion

info <- matrix(NA, nrow = 2, ncol = 2)
info[2,1] <- my_aic(likelihood(full_sample[["dccoutput"]][["fit"]]),13)
info[2,2] <- my_aic((likelihood(full_sample[["dccoutput_tree_2"]][["fit"]])+likelihood(full_sample[["dccoutput_tree_1"]][["fit"]])),13)
info[1,1] <- likelihood(full_sample[["dccoutput"]][["fit"]])
info[1,2] <- likelihood(full_sample[["dccoutput_tree_1"]][["fit"]])+likelihood(full_sample[["dccoutput_tree_2"]][["fit"]])



rownames(info) <- c("LogLik", "AIC")
colnames(info) <- c("DCC-GARCH", "DCC-Tree-GARCH")
xtable(info,digits = 5)