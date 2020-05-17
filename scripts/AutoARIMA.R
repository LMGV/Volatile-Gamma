library(forecast)
library(tseries)
setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub/data")

data <- read.table("data.csv", sep = ",")

#### ACF and PACF ####
acf(data$oil, lag.max = 30, plot = TRUE)
acf(data$rub, lag.max = 30, plot = TRUE)

data2 <- data^2
acf(data2$oil, lag.max = 30, plot = TRUE)
acf(data2$rub, lag.max = 30, plot = TRUE)

#### Ljung-Box test ####

# H0: independence
Box.test(data$oil, lag = 2, type = "Ljung-Box")#$p.value
Box.test(data$rub, lag = 2, type = "Ljung-Box")#$p.value
#auto.arima(data$oil, stationary = TRUE, ic = "bic")
#auto.arima(data$rub, stationary = TRUE, ic = "bic")

#### Dickey-Fuller test ####

# H0: unit root
#adf.test(x = data$oil, k = 0)#$p.value

oil.t   <- data$oil[-1]
oil.t_1 <- data$oil[-length(data$oil)]
x.oil   <- matrix(c(rep(1, length(oil.t_1)), oil.t_1), nrow = length(oil.t_1), ncol = 2) # matrix of a constant and lagged dep var

rub.t   <- data$rub[-1]
rub.t_1 <- data$rub[-length(data$rub)]
x.rub   <- matrix(c(rep(1, length(rub.t_1)), rub.t_1), nrow = length(rub.t_1), ncol = 2) # matrix of a constant and lagged dep var


######### fit regression y_t = \alpha0 + \alpha1*y_{t-1} + \epsilon_t

# this function performs a DF test with a constant; returns 1 if differencing is required and 0 otherwise
df.test <- function(y, X, cl){ # cl can be one of 0.01, 0.025, 0.05, 0.1
  if( !(cl %in% c(0.01, 0.025, 0.05, 0.1)) ) stop('cl must be one of 0.01, 0.025, 0.05, 0.1')
  betas.hat   <- solve(crossprod(X)) %*% t(X) %*% y
  sigma.2.hat <- crossprod(y - (X %*% betas.hat)) / (nrow(X) - ncol(X))
  varcovm     <- as.vector(sigma.2.hat) * solve(crossprod(X))
  s.e         <- sqrt(diag(varcovm))
  d.f.stat    <- (betas.hat[2]-1)/s.e[2]
  conf.levels <- c(0.01, 0.025, 0.05, 0.1)
  crit.values <- c(-3.434, -3.120, -2.863, -2.568)
  
  if (d.f.stat < crit.values[match(cl, conf.levels)]){
    d <- 0
  } else {
    d <- 1
  }
  
  return(d)
}

d.oil <- df.test(oil.t, x.oil, 0.05)
d.rub <- df.test(rub.t, x.rub, 0.05)


#### ARIMA model fit ####

# this function computes sum of squared errors of a given ARMA model
arma.errors <- function(parameters, y, p, q){ # parameters should be passed as a vector of c, AR, MA, sigma^2
  if( any(p < 0 | q < 0 | p%%1!=0 | q%%1!=0) ) stop('p and q must be non-negative integers')
  epsilon    <- rep(0, length(y))
  epsilon[1] <- y[1] - parameters[1]
  
  # to allow for varying p and q
  
  ar.term <- function(parameters, y, p, o){
    if (p == 0){return(0)}
    ar.vect <- rep(0, p)
    for (j in 1:p){
      ar.vect[j] <- parameters[1+j]*y[o-j]
    }
    return(sum(ar.vect))
  }
  
  ma.term <- function(parameters, e, q, o, p){
    if (q == 0){return(0)}
    ma.vect <- rep(0, q)
    for (k in 1:q){
      ma.vect[k] <- parameters[1+p+k]*e[o-k]
    }
    return(sum(ma.vect))
  }
  
  if (p > 1){
    if (q > 1){
      if (p < q){
        
        for (i in (2:p)){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, i-1, i) - ma.term(parameters, epsilon, i-1, i, p)
        }
        for (i in ((p+1):q)){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, i-1, i, p)
        }
        for (i in ((q+1):length(y))){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
        }
        
      }
      
      if (p > q){
        
        for (i in (2:q)){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, i-1, i) - ma.term(parameters, epsilon, i-1, i, p)
        }
        for (i in ((q+1):p)){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, i-1, i) - ma.term(parameters, epsilon, q, i, p)
        }
        for (i in ((p+1):length(y))){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
        }
        
      }
      
      if (p == q){
        
        for (i in (2:q)){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, i-1, i) - ma.term(parameters, epsilon, i-1, i, p)
        }
        for (i in ((q+1):length(y))){
          epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
        }
        
      }
      
      
    }
    
    else{ # q = 1 or q = 0, p > 1
      
      for (i in (2:p)){
        epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, i-1, i) - ma.term(parameters, epsilon, q, i, p)
      }
      for (i in ((p+1):length(y))){
        epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
      }
      
    }
  }
  else{ # p = 1 or p = 0
    if (q > 1){
      
      for (i in (2:q)){
        epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, i-1, i, p)
      }
      for (i in ((q+1):length(y))){
        epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
      }
      
    }
    
    else{ # p = 1, q = 1, p = 0, q = 0
      
      for (i in (2:length(y))){
        epsilon[i] <- y[i] - parameters[1] - ar.term(parameters, y, p, i) - ma.term(parameters, epsilon, q, i, p)
      }
      
    }
  }
  
  return(sum(epsilon^2))
  
}

# this function computes log-likelihood value of a given ARMA model
log.lik <- function(parameters, y, p, q) {
  e  <- arma.errors(parameters, y, p, q)
  ll <- -(length(y)/2)*log(2*pi) - (length(y)/2)*log(parameters[length(parameters)]) - (1/(2*parameters[length(parameters)]))*e
  return(-1 * ll)
}

# this function fits an ARIMA model and returns estimated coefficients, value of LL function and value of BIC
fit.ARIMA <- function(y, p, d, q){ # data vector, AR order, differencing, MA order
  if (d == 1){
    y = diff(y)
  }
  
  arma.coeffs  <- optim(par = c(rep(0, p+q+1), 1), fn = log.lik, y = y, p = p, q = q)$par # minimize the negative log-likelihood
  loglik.value <- -1*log.lik(arma.coeffs, y, p, q)                                        # get the log-likelihood value
  bic          <- -2*loglik.value + (p + q + 1)*log(length(y))                            # compute the BIC
  
  returns.list  <- list("Coefficients" = arma.coeffs, "LogLikelihoodvalue" = loglik.value, "BIC" = bic)
  
  return(returns.list)
  
}

# this function returns the best ARIMA model according to the BIC criterion
best.ARIMA <- function(y, p_grid, d, q_grid){
  model.bic    <- 0
  
  for (pv in p_grid){
    for (qv in q_grid){
      print(paste("Trying model:", pv, d, qv))
      model <- fit.ARIMA(y, pv, d, qv)
      print(paste("BIC:", model$BIC))
      if (model$BIC < model.bic){
        model.bic    <- model$BIC
        model.coeffs <- model$Coefficients
        model.loglik <- model$LogLikelihoodvalue
        model.p      <- pv
        model.q      <- qv
      }
    }
  }
  
  best.model <- list("Coefficients" = model.coeffs, "LogLikelihoodvalue" = model.loglik, "BIC" = model.bic, "p" = model.p, 
                     "q" = model.q, "d" = d)
  return(best.model)
}

model.oil <- best.ARIMA(oil.t, c(0, 1, 2), d.oil, c(0, 1, 2))
model.rub <- best.ARIMA(rub.t, c(0, 1, 2), d.rub, c(0, 1, 2))

