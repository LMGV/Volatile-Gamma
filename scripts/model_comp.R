QLIKE <- function(returns, estimates){
  a=0
  for (i in 1:length(mod[,1,1])){
    a=a-log(dmvnorm(returns[i,4:7],mean=colMeans(returns[i,4:7]),sigma=estimates[i,4:7,4:7]))
  }
  a
}
####Funciton for real.cov 
real_cov <- function(returns){#only works for bivariate models
  real_cov <- array(NA,dim = c(nrow(returns),2,2))
  for (i in 1:nrow(returns)) {
    real_cov[i,2,1] <- (mean(returns[,1])-returns[i,1])*(mean(returns[,2])-returns[i,2])
    real_cov[i,1,2] <- real_cov[i,2,1]
    real_cov[i,1,1] <- (mean(returns[,1])-returns[i,1])*(mean(returns[,1])-returns[i,1])
    real_cov[i,2,2] <- (mean(returns[,2])-returns[i,2])*(mean(returns[,2])-returns[i,2])
    }
  real_cov
}


MEA_MSE <- function(real_cov, estimates){
  a=matrix(,4,4)  #MAE/MSE
  b=matrix(,4,4)
  
  for (i in 4:7)
  {
    for (j in 4:7)
    {   
      a[i-3,j-3]=mean(abs(real_cov[,i,j]-estimates[,i,j]))
      a[j-3,i-3]=a[i-3,j-3]
      b[i-3,j-3]=mean((real_cov[,i,j]-estimates[,i,j])^2)
      b[j-3,i-3]=b[i-3,j-3]
    }
  }
  
  MAE <- mean(a)
  MSE <- mean(b)
  output <- rbind(MAE,MSE)
}

DMW <- function(returns, estimates_1, estimates_2, real_cov, loss_function){ #1= QLIKE, 2= MAE, 3=MSE
  n=length(estimates_1[,1,1])
  
  perf=rep(0,length(estimates_1[,1,1])) 
  
  for (i in 1:length(estimates_1[,1,1])){
    
    if (loss_function==1) {
      perf[i]= length(estimates_1[,1,1])*(-log(dmvnorm(returns[i,4:7],mean=colMeans(returns[i,4:7]),sigma=estimates_1[i,4:7,4:7]))+log(dmvnorm(returns[i,4:7],mean=colMeans(returns[,4:7]),sigma=estimates_2[i,4:7,4:7]))) #QLIKE
    } else if (loss_function==2) {
      
      perf[i]= mean(abs(real_cov[i,,]-us.stock.data$ccc.estimates[i,4:7,4:7])) - mean(abs(real_cov[i,4:7,4:7]-us.stock.data$dcc.estimates[i,4:7,4:7])) #MAE
    
      } else if (loss_function==3) {
    perf[i]= mean(abs(real_cov[i,,]-us.stock.data$ccc.estimates[i,4:7,4:7])^2) - mean(abs(real_cov[i,4:7,4:7]-us.stock.data$dcc.estimates[i,4:7,4:7])^2) #MSE
    } else {}
  }
  DMW_test <- sqrt(n)*mean(perf)/sqrt(spectrum(perf)$spec[1]) 
  DMW_p <- 1-pnorm(sqrt(n)*mean(perf)/sqrt(spectrum(perf)$spec[1]))
  DMW <- cbind(DMW_test,DMW_p) 
}
  
  
model_comparison <- function(estimates_1, estimates_2, returns, loss_function){
  QLIKE_1 <- QLIKE(returns, estimates_1)
  QLIKE_2 <- QLIKE(returns, estimates_2)
  loss_function <- loss_function
  real_cov <- NA
  #real_cov <- real_cov(returns)
  #MAEMSE_1 <- MEA_MSE(real_cov, estimates_1)
  #MAEMSE_2 <- MEA_MSE(real_cov, estimates_2)
  DMW <- DMW(returns, estimates_1, estimates_2, real_cov, loss_function)
  comparison <- matrix(data = NA,nrow = 3,ncol = 2)
  #comparison[1:2,1] <- MAEMSE_1
  #comparison[1:2,2] <- MAEMSE_2
  comparison[3,1] <- QLIKE_1
  comparison[3,2] <- QLIKE_2
  out_put <- list("Loss_function"= comparison, "DMW"= DMW)
  out_put
}

estimates_1=us.stock.data$ccc.estimates ##DCC model estimates length(ts) * 2 * 2
estimates_2=us.stock.data$dcc.estimates

returns <- us.stock.data$returns

output <- model_comparison(estimates_1, estimates_2, returns, loss_function = 1)
