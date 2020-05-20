dccforecast(dccfit, n.ahead = 1, n.roll = nrow(returns))


dccforecast <- function(dccfit, returns){
  # get all the  D_ts from eric
  D_t <-  matrix(data = 0, nrow = 2, ncol = 2) #Diagonal Matrix with the sigma_t_i
  #get all the vts
  vt <- inv(D_t) * (r_t - mu)# (r_t - mu) row vector
  
  #Constants 
  dcca1 <- dccfit@model[["pars"]][1,1]
  dccb1 <- dccfit@model[["pars"]][2,1]
  R <- 1/nrow(returns) sum(v_t1 %*% t(v_ti))
  #Changeing over t
  Q_t <- rep(list(NA),nrow(returns)
  for (i in 2:nrow(returns)) {
  Q_t_1 <- dccfit@mfit[["Q"]][[i-1]] #get erics variance predictions in there
  
  v_t <- inv(D_t[i-1]) * (r_t[i-1] - mu)# (r_t - mu) row vector
  
  
  Q_t[[i]] <- R + dcca1(v_t %*% t(v_t) - R) + dccb1(Q_t_1 - R) 
  }
}

