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

Box.test(data$oil, lag = 100, type = "Ljung-Box")#$p.value
auto.arima(data$oil, stationary = TRUE, ic = "bic")
auto.arima(data$rub, stationary = TRUE, ic = "bic")

#### Dickey-Fuller test ####

# H0: unit root
#adf.test(x = data$oil, k = 0)#$p.value

oil.t   <- data$oil[-1]
oil.t_1 <- data$oil[-length(data$oil)]
x       <- matrix(c(rep(1, length(oil.t_1)), oil.t_1), nrow = length(oil.t_1), ncol = 2) # matrix of a constant and lagged dep var

# fit regression y_t = \alpha0 + \alpha1*y_{t-1} + \epsilon_t
betas.hat   <- solve(crossprod(x)) %*% t(x) %*% oil.t
sigma.2.hat <- crossprod(oil.t - (x %*% betas.hat)) / (nrow(x) - ncol(x))
varcovm     <- as.vector(sigma.2.hat) * solve(crossprod(x))
s.e         <- sqrt(diag(varcovm))
d.f.stat    <- (betas.hat[2]-1)/s.e[2]
crit.values <- c(-3.434, -3.120, -2.863, -2.568) # Dickey-Fuller critical values for 1%, 2.5%, 5% and 10% respectively
if (d.f.stat < crit.values[3]){                  # choose 5% significance level
  d <- 0
} else {
  d <- 1
}


# diff(y, differences = 2) - second difference


#### ARMA model fit ####

# p = 1, q = 1
# 1. Estimate model parameters (OLS/MLE)
# 2. Get Likelihood of the model
# 3. Compute BIC or AICc

# 
# regression: y_t = c + \epsilon_t + \phi_1*y_{t-1} + \theta_1\eta_{t-1}


arma.errors <- function(parameters, y){ # parameters should be passed as a vector of c, AR, MA
  epsilon <- rep(0, length(y))
  epsilon[1] <- y[1] - parameters[1]
  
  for (i in (2:length(y))){
    epsilon[i] <- y[i] - parameters[1] - parameters[2]*y[i-1] - parameters[3]*epsilon[i-1]
  }
  
  return(sum(epsilon^2))
  
}

#### some something (delete later) ####
#arma.errors <- function(parameters){ # parameters should be passed as a vector of c, AR, MA
#  epsilon <- rep(0, length(oil.t))
#  epsilon[1] <- oil.t[1] - parameters[1]
#  
#  for (i in (2:length(oil.t))){
#    epsilon[i] <- oil.t[i] - parameters[1] - parameters[2]*oil.t[i-1] - parameters[3]*epsilon[i-1]
#  }
#  
#  return(sum(epsilon^2))
#  
#}

#arma.coeffs <- nlm(f = arma.errors, p = c(0, 0, 0), y = oil.t)$gradient
#arma.coeffs <- nlminb(objective = arma.errors, start = c(0, 0, 0))$par
#arma.coeffs <- optim(par = c(0,0,0), fn = arma.errors, y = oil.t)$par

#epsilon_mle <- rep(0, length(oil.t))
#epsilon_mle[1] <- oil.t[1] - arma.coeffs[1]
#for (i in (2:length(oil.t))){
#  epsilon_mle[i] <- oil.t[i] - arma.coeffs[1] - arma.coeffs[2]*oil.t[i-1] - arma.coeffs[3]*epsilon_mle[i-1]
#}

#log(sum(epsilon_mle^2))
#### ####

# fit a log-likelihood function
log.lik <- function(parameters, y) {
  e  <- arma.errors(parameters, y)
  ll <- -(length(y)/2)*log(2*pi) - (length(y)/2)*log(parameters[4]) - (1/(2*parameters[4]))*e
  return(-1 * ll)
}

arma.coeffs  <- optim(par = c(0,0,0,1), fn = log.lik, y = oil.t)$par # minimize the negative log-likelihood
loglik.value <- -1*log.lik(arma.coeffs, oil.t)                       # get the log-likelihood value
p = 1
q = 1
bic          <- -2*loglik.value + ((p + q + 1))*log(length(oil.t))   # compute the BIC



