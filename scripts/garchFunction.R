# garch function

garchEstimation = function(theta, returns, ar, ma, threshhold,data_threshhold,type, distribution) { 
  {
    # Inputs
      # theta: AR-coefs, MA-coefs, Treshhold, df of t-distribution
      # returns: univariate xts vector
      # ar: order AR process for squared returns
      # ma: order MA process for squared returns
      # threshhold: T/F. if T, then # cols of data_treshhold is number of threshhold parameters. if F then inactive
      # data_treshhold: time series object. only evaluated if threshhold = T
      # type: GARCH-model type.
      # distribution: normal / t
    
    # Outputs
      # likelihood of sample given the model specification
    
    # conditions
       # check input is xts
      if (is.xts(returns)==F) {
        print("Error: No Time Series Object as input to garchEstimation")
      }
      # distribution
      if (distribution %in% c("normal", "norm", "t") ==F ){
        print("Error: non-supported distribution type in garchEstimation")
      } 
      # number of parameters to be estimated equals dim of theta
      if (length(theta) != (1+1 + ar+ma+as.numeric(threshhold) + as.numeric(distribution=="t"))) { # number of parms: mean return+ constant + ar +ma + threshhold(if active) + degrees of freedom t (if active)
        print("Error: Number of input parameters does not match length of parameter vector in garchEstimation")
      } 


    n=length(returns)
    x.start= mean(returns)
    sigmasq.start= var(returns)
    
    data=c(x.start,returns)
    sigmasq= rep(0,n+1)
    sigmasq[1]=sigmasq.start
    
    my.sigma=c(sqrt(sigmasq[1]),rep(0,n))
    
    mean_ret=rep(theta[1],n+1) # constant conditional mean
    
    for (i in 2:(n+1))
    {
      ar_part = 
      ma_part  =
      sigmasq  =theta[2]^2 # force constant to be positive
      sigmasq= theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]^2*sigmasq[i-1] #GARCH(1,1)
      sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]^2*sigmasq[i-1] #GARCH(1,1)
      
      #sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])^2 + theta[4]*(data[i-1]-mean_ret[i-1])^2*((data[i-1]-mean_ret[i-1])<=0)+theta[5]^2*sigmasq[i-1] #GJR-GARCH(1,1)
      
      #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-mean_ret[i-1])*((data[i-1]-mean_ret[i-1])&gt;0)-theta[4]^2*(data[i-1]-mean_ret[i-1])*((data[i-1]-mean_ret[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
      
      #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-mean_ret[i-1]))-theta[4]*(data[i-1]-mean_ret[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
      #sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-mean_ret[i-1]))-theta[4]*(data[i-1]-mean_ret[i-1]))^2+theta[5]^2*sigmasq[i-1] #PGARCH(1,1) with d=2
    }
    
    # normdistrib, GARCH 1/1
    1/2*sum(log(sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]))))
    #tdistrib
    # 1/2*sum(log(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
  }
}

