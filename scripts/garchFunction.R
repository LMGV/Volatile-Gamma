# garch function

garchEstimation = function(theta, returns, ar, ma, threshhold,th_value,data_threshhold, distribution) { 
  {
    # Inputs
      # theta: AR-coefs, MA-coefs, Treshhold, df of t-distribution
      # returns: univariate return vector, double
      # ar: order AR process for squared returns
      # ma: order MA process for squared returns
      # threshhold: T/F. if T, then # cols of data_treshhold is number of threshhold parameters. if F then inactive
      # data_treshhold: only evaluated if threshhold = T. ATM not used
      # distribution: normal / t
    
    # Outputs
      # likelihood of sample given the model specification
      # !parameters square root! (when using an optimizer)
    
    # conditions
      # distribution
      if (distribution %in% c("normal", "t") ==F ){
        print("Error: non-supported distribution type in garchEstimation")
      } 
      # number of parameters to be estimated equals dim of theta
      if (length(theta) != (1+1 + ma+ ar+ as.numeric(threshhold) + as.numeric(distribution=="t"))) { # number of parms: mean return+ constant + ar +ma + threshhold_parameter (if active) + degrees of freedom t (if active)
        print("Error: Number of input parameters does not match length of parameter vector in garchEstimation")
      } 
    
    
      # assign coefficients:
      mu_coef = theta[1]
      constant_coef = theta[2]
      ma_coef = theta[(1:ma)+2]
      ar_coef = theta[((ma+1):(ma+ar))+2]
      th_coef = ifelse(threshhold==T, theta[(ar+ma+1)+2],NA)
      df_t_coef = ifelse(distribution=="t", theta[(ar+ma+1+as.numeric(threshhold))+2],NA)
    
      # number of timesteps, starting values
      max_lags = max(ar, ma)
      n=length(returns)
      x.start= mean(returns)
      sigmasq.start= var(returns)
      
      # define returns
      data=c(x.start,returns)
      
      # initialize variance as sample variance
      sigmasq= rep(0,n+1)
      sigmasq[1]=sigmasq.start # initialize with unconditional variance
      
      my.sigma=c(sqrt(sigmasq[1]),rep(0,n))
      
      mean_ret=rep(mu_coef,n+1) # constant conditional mean
      # mean_ret[1] = 0 
      ar_part = rep(0, n+1) # initialize AR
      ma_part = rep(0, n+1) # initialize MA
      th_active_iteration =  rep(0, n+1) # initialize th_part (0 if not active)
      
      # calculate ma part. loop for all ma parts
      for (k in 1:length(ma_coef)) {
        ma_part[(1+k):(n+1)] = ma_part[(1+k):(n+1)] + ma_coef[k]^2*(data[1:(n+1-k)]-mean_ret[1:(n+1-k)])^2 # calc MA part
      }

      # threshhold: th_value as input parameter is not optimized. Parm "active" when epsilon < th_value
      if (threshhold==T)  {
        th_active_iteration[2:(n+1)] = ((data[1:n]-mean_ret[1:n])<=th_value) * ((data[1:n]-mean_ret[1:n])^2)*th_coef
      }

      # loop for time: calc AR and conditional variance
      for (i in (max_lags+1):(n+1))
      {
        # calculate AR parts. loop for all AR parts
        for (j in 1:length(ar_coef)) {
          ar_part[i] = ar_part[i] + ar_coef[j]^2*sigmasq[i-j]  # calc AR part
        }
        # ar_part[i]  = theta[3]^2*sigmasq[i-1] # iteratively add ar part
        sigmasq[i]  = constant_coef^2 + ar_part[i] + ma_part[i] +th_active_iteration[i] # conditional variance. force constant to be positive
      }
    
    # calc log likelihood of model
      if (distribution=="normal"){
        log_liklihood = (n+1-max_lags)*log(sqrt(2*pi))+sum(0.5*((data[(max_lags+1):(n+1)]-mean_ret[(max_lags+1):(n+1)])^2)/sigmasq[(max_lags+1):(n+1)]) + sum(0.5*log(sigmasq[(max_lags+1):(n+1)]))
        return(log_liklihood)
      } else  if (distribution=="t"){
        log_liklihood = -(n+1-max_lags)*log(gamma((df_t_coef+1)/2)/(gamma(df_t_coef/2)*sqrt(pi*(df_t_coef-2)))) + (df_t_coef+1)/2*sum(log(1+((data[(max_lags+1):(n+1)]-mean_ret[(max_lags+1):(n+1)])^2)/((df_t_coef-2)*sigmasq[(max_lags+1):(n+1)]))) + 0.5*sum(log(sigmasq[(max_lags+1):(n+1)]))
        return(log_liklihood)
      }
      
    # audrino implementation likelihood to compare. Delivers same results
      # normal
           # 1/2*sum(log(sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)])))) # audrino code
      #tdistrib 
           # 1/2*sum(log(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10) # audrino code
  }
}


my.loglike.t=function(theta) #Estimate an asymmetric GARCH(1,1) model with Student's t innovations
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
  
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start
  
  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,n))
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }
  
  for (i in 2:(n+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)
    
    # my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    # my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }
  
  #my.sigmasq=my.sigma^2cd
  #my.sigmasq=exp(log.sigmasq)
  
  # normdistrib, GARCH 1/1
  1/2*sum(log(my.sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]))))
  
  #tdistrib
  #1/2*sum(log(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
  # return(-n*log(gamma((theta[6]+1)/2)/(gamma(theta[6]/2)*sqrt(pi*(theta[6]-2)))) + (theta[6]+1)/2*sum(log(1+((data[2:(n+1)]-my.mean[2:(n+1)])^2)/((theta[6]-2)*my.sigmasq[2:(n+1)]))) + 0.5*sum(log(my.sigmasq[2:(n+1)])))
}
