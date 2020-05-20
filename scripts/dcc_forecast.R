library(fBasics)


dcc_forecast <- function(dccfit, pred_results_1, pred_results_2){
  D <- array(0,dim = c(length(pred_results_1),2,2))
  for (i in 1:length(pred_results_1)) {
    D[i,1,1] <- sqrt(pred_results_1[i])
    D[i,2,2] <- sqrt(pred_results_2[i])
  }
  V <- array(0,dim = c(length(pred_results_1),2,2))
  for (i in 1:length(pred_results_1)) {
    V[i,,] <- inv(D[i,,])#%*%(returns[i,] - colMeans(returns))# add retruns tomorrow
  }

  dcca1 <- dccfit@model[["pars"]][1,1]
  dccb1 <- dccfit@model[["pars"]][2,1]
  
  VtV <- array(0,dim = c(length(pred_results_1),2,2))
  for (i in 1:length(pred_results_1)) {
    VtV[i,,] <- V[i,,] %*% t(V[i,,])
  }
  
  R <- sum(VtV)/length(pred_results_1)
  
  Q_t <- array(0,dim = c(length(pred_results_1),2,2))
  for (i in 2:length(pred_results_1)) {
    Q_t[i,,] <- R + dcca1*(VtV[i,,] - R) + dccb1*(Q_t[i-1,,] - R) #What is the initial Q_t
  }
  Q_t
}

pred_results <- readRDS("C:/Users/johan/Documents/GitHub/Volatile-Gamma/output/univariateModels/in_sample_pred_result.rds")
pred_rub <- pred_results$rub$variance_predict
pred_oil <- pred_results$oil$variance_predict

estimates <- dcc_forecast(dccfit1,pred_rub,pred_oil)
