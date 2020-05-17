# garch function

garchEstimation = function(theta, returns, ar, ma, threshhold,th_value,data_threshhold,type, distribution) { 
  {
    # Inputs
      # theta: AR-coefs, MA-coefs, Treshhold, df of t-distribution
      # returns: univariate return vector, double
      # ar: order AR process for squared returns
      # ma: order MA process for squared returns
      # threshhold: T/F. if T, then # cols of data_treshhold is number of threshhold parameters. if F then inactive
      # data_treshhold: only evaluated if threshhold = T
      # type: GARCH-model type.
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
      if (length(theta) != (1+1 + ar+ma+as.numeric(threshhold) + as.numeric(distribution=="t"))) { # number of parms: mean return+ constant + ar +ma + threshhold_parameter (if active) + degrees of freedom t (if active)
        print("Error: Number of input parameters does not match length of parameter vector in garchEstimation")
      } 

      # number of timesteps, starting values
      max_lags = max(ar, ma)
      n=length(returns)
      x.start= mean(returns)
      sigmasq.start= var(returns)
      
      # define returns
      data=c(x.start,returns)
      
      # initialize variance as sample variance
      sigmasq= rep(0,n+1)
      sigmasq[1]=sigmasq.start
      
      my.sigma=c(sqrt(sigmasq[1]),rep(0,n))
      
      mean_ret=rep(theta[1],n+1) # constant conditional mean
      # mean_ret[1] = 0 
      ar_part = rep(0, n+1) # initialize AR
      ma_part = rep(0, n+1) # initialize MA
      ma_part[2:(n+1)] = theta[4]^2*(data[1:n]-mean_ret[1:n])^2 # calc MA part
      th_active_iteration =  rep(0, n+1) # initialize th_part (0 if not active)
      
      # threshhold: th_value as input parameter is not optimized. Parm "active" when epsilon < th_value
      if (threshhold==T)  {
        th_active_iteration[2:(n+1)] = ((data[1:n]-mean_ret[1:n])<=th_value) * ((data[1:n]-mean_ret[1:n])^2)*theta[5]
      }

      for (i in 2:(n+1))
      {
        ar_part[i]  = theta[3]^2*sigmasq[i-1] # iteratively add ar part
   
        sigmasq[i]  =theta[2]^2 + ar_part[i] + ma_part[i] +th_active_iteration[i]# force constant to be positive
        # sigmasq[i]= theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]^2*sigmasq[i-1] #GARCH(1,1)
        # sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]^2*sigmasq[i-1] #GARCH(1,1)
        
        # sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]*(data[i-1]-mean_ret[i-1])^2*((data[i-1]-mean_ret[i-1])<=0)+theta[5]^2*sigmasq[i-1] #GJR-GARCH(1,1)
        
        #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])*((data[i-1]-mean_ret[i-1])&gt;0)-theta[4]^2*(data[i-1]-mean_ret[i-1])*((data[i-1]-mean_ret[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
        
        #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-mean_ret[i-1]))-theta[4]*(data[i-1]-mean_ret[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
        #sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-mean_ret[i-1]))-theta[4]*(data[i-1]-mean_ret[i-1]))^2+theta[5]^2*sigmasq[i-1] #PGARCH(1,1) with d=2
      }
    
    # calc log likelihood of model
      
    # normdistrib, GARCH 1/1 (same for GJR)
    1/2*sum(log(sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]))))
    #tdistrib
    # 1/2*sum(log(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
  }
}

