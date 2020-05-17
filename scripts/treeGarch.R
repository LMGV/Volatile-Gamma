# Tree Garch - Test

#### Pre
# librarys
  library(readxl)    
  library(readr)   
  library(data.table)
  library(tidyverse)  
  library(stats)
  library(ggfortify)  
  library(zoo)        
  library(xtable)     
  library(lmtest)     
  library(forecast) 
  library(tseries)
  library(xts)
  library(quantmod)
  library(fGarch)
  library(rugarch)
  library(rmgarch)
  library(psych)
  library(MASS)
  library(fitdistrplus)
  
  filter <- dplyr::filter
  select <- dplyr::select
  setwd("~/GitHub/Volatile-Gamma") # setwd
source("scripts/functions.R") # functions
source("scripts/garchFunction.R") # functions

# parameters
quantile_outliers = 0.001 # cut 0.1% of right and left tail


####Data Import####
  ts_r = read.table('data/data.csv', sep = ',')
  ts_r = xts(ts_r, order.by = as.Date(rownames(ts_r)))
  
  
# descriptives / arima. replace later ----
          # Remove Outliers via empirical quantiles
          ts_r <- ts_r[(ts_r$rub > quantile(ts_r$rub,quantile_outliers)) & (ts_r$rub < quantile(ts_r$rub,1-quantile_outliers)) & (ts_r$oil > quantile(ts_r$oil,quantile_outliers))& (ts_r$oil < quantile(ts_r$oil,1-quantile_outliers))]
          
          
          ####Mean####
          print("Mean of Rubel log returns")
          mean(ts_r$rub)
          
          print("Mean of oil log returns")
          mean(ts_r$oil)
          
          ####Summary Statitics####
          summary_oil <- describe(ts_r$oil)
          summary_rub <- describe(ts_r$rub)
          
          ####Test for Normalily QQ plot####
          qqnorm(ts_r$rub)
          qqline(ts_r$rub)
          
          qqnorm(ts_r$oil)
          qqline(ts_r$oil)
          
          jarque.bera.test(ts_r$rub)
          jarque.bera.test(ts_r$oil)
          
          ####Fit####
          # CHECK which method returns plausible results. vary outlier quantile and see if t fits better then
          fitdist(as.vector(ts_r$rub),"t", method = "mge", start = list(df=2))
          
          ####ACF r####
          acf(ts_r$oil, lag.max = 30, plot = TRUE)
          acf(ts_r$rub, lag.max = 30, plot = TRUE)
          
          pacf(ts_r$rub, lag.max = 30, plot = TRUE)
          pacf(ts_r$oil, lag.max = 30, plot = TRUE)
          
          ####ACF r^2####
          ts_r2 <- ts_r^2
          acf(ts_r2$oil, lag.max = 30, plot = TRUE)
          acf(ts_r2$rub, lag.max = 30, plot = TRUE)
          
          ####ARIMA####
          arma_oil <- auto.arima(ts_r$oil)
          arma_rub <- auto.arima(ts_r$rub)
          
          oil_ar <- arma_oil[["arma"]][1]
          oil_ma <- arma_oil[["arma"]][2]
          
          rub_ar <- arma_oil[["arma"]][1]
          rub_ma <- arma_oil[["arma"]][2]

          

# Garch Function ----
          source("scripts/garchFunction.R") # functions
          # inputs fct
          returns=ts_r$oil
          ar = 1
          ma = 1
          threshhold = F
          th_value  = 0 # not optimized within fct
          data_threshhold = 0
          type = "GARCH"
          distribution ="norm"
          
          # my function 
          start_parms=c(rep(0.5,4),0.1)
          opt_parms= nlm(garchEstimation,start_parms,
                         returns = returns,  ar = ar, ma = ma,
                         threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                         type=type, distribution=distribution,
                         print.level=2,iterlim=1000, check.analyticals=1)
          
          names_coefs = c("exp_ret","constant", paste0("AR",c(1:ar)), paste0("MA",c(1:ar)), "threshhold_parm")
          names_coefs
          c(opt_parms$estimate[1], opt_parms$estimate[2:(2+ar+ma)]^2,opt_parms$estimate[(3+ar+ma)])           # make sure to revert squares in parms
            
          # news impact curve
            
          garchEstimation(theta=start_parms, returns, ar, ma, threshhold,th_value,data_threshhold,type, distribution)
          my.loglike.t(start_parms)
          
          # Audrino fct
          par.start=c(rep(0.5,4),0.1)
          my.optpar= nlm(my.loglike.t,par.start,iterlim=1000,print.level=1)
          my.optpar$estimate
          
          names_coefs = c("exp_ret","constant", "MA1",  "threshhold_parm", "AR1")
          names_coefs
          c(my.optpar$estimate[1], my.optpar$estimate[2]^2,my.optpar$estimate[3]^2,my.optpar$estimate[4] ,my.optpar$estimate[5]^2 )          # make sure to revert squares in parms
          
          
          m1=garchFit(returns~garch(1,1),data=returns,trace=F)
          summary(m1)
          
          # ug_spec= ugarchspec(variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(0, ), include.mean = TRUE), distribution.model="norm")
          # ug_fit= ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
          # ug_fit
          
          
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
    
              #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
              
              #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
              #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
            }
            
            #my.sigmasq=my.sigma^2cd
            #my.sigmasq=exp(log.sigmasq)
            
            # normdistrib, GARCH 1/1
            1/2*sum(log(my.sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]))))
            #tdistrib
            # 1/2*sum(log(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
          }
          
          ####Univariate Garch Opt function ####
          opt_garch<- function(ar,ma,returns) {
            AIC <- matrix(data=NA,nrow=2,ncol=2)
            for (j in 1:2) {
              for (i in 1:2) {
                ug_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="sstd") #,submodel="TGARCH"
                ugfit = ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
                info <- infocriteria(ugfit)
                AIC[j,i] <- info[2,1]
              }
            }
            opt_garch <- which(AIC == min(AIC,na.rm=TRUE), arr.ind=TRUE)
            ug_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(opt_garch[1],opt_garch[2])), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="sstd") # ,submodel="TGARCH"
            ug_fit = ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
            output <- list("Specs"=ug_spec,"fit"=ug_fit)
          }
          
          ####Uni-Garch Specs Oil and RUB####
          output_oil <- opt_garch(oil_ar,oil_ma, ts_r$oil)
          output_rub <- opt_garch(rub_ar,rub_ma, ts_r$rub)
          
          #Garch summary Oil
          output_oil[["fit"]]
          
          #Garch summary Oil
          output_rub[["fit"]]
          
          ####Bivar Garch####
          cor(ts_r$rub,ts_r$oil)
          var(ts_r)
          var(ts_r^2)
          cor(ts_r$rub^2,ts_r$oil^2)
          
