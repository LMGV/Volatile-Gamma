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
          
          ####ARIMA####
          arma_oil <- auto.arima(ts_r$oil)
          arma_rub <- auto.arima(ts_r$rub)
          
          oil_ar <- arma_oil[["arma"]][1]
          oil_ma <- arma_oil[["arma"]][2]
          
          rub_ar <- arma_oil[["arma"]][1]
          rub_ma <- arma_oil[["arma"]][2]

          

# Garch Function ----
    # to see number of lags: just do 
          
    # ACF / PACF r^2 
    #need to adjust confidence bounds, portm-test
      ts_r2 <- ts_r^2
      acf(ts_r2$oil, lag.max = 30, plot = TRUE)
      acf(ts_r2$rub, lag.max = 30, plot = TRUE)  
      pacf(ts_r2$oil, lag.max = 30, plot = TRUE)
      pacf(ts_r2$rub, lag.max = 30, plot = TRUE)  
      # suggests high order AR lags
          
          
    # test for asymmetries
      # autocorrelation
      sampleAutocorrelation(ts_r$oil, "oil", 0.05)
      sampleAutocorrelation(ts_r$rub, "RUBUSD", 0.05)

          
     # sign bias tests   
     # TODO
      a=(sp500<0) #Sign Bias
      a=c()
      for (i in 1:length(ibm))
      {
        a=c(a,100*min(sp500[i],0)) #Negative Size Bias
        #a=c(a,100*max(sp500[i],0)) #Positive Size Bias
      }

      h=5
      b=lm(100*sp500[(h+1):length(ibm)]^2~a[1:(length(ibm)-h)])
      summary(b)
          
          
     
   # tree garch test ----
            returns = ts_r$oil

              # define possible split variables
                # past returns, epsilons, variances of own process and other processes
                 means = colMeans(ts_r)
                 epsilon = sweep(ts_r,2,means)
                 epsilon_sq = epsilon^2
                 colnames(epsilon) = paste0(colnames(ts_r),"_epsilon")
                 colnames(epsilon_sq) = paste0(colnames(ts_r),"_epsilon_sq")
                 
                 
                 base_split_variables = as.xts(cbind(ts_r, epsilon, epsilon_sq))
                 
                 # get 2 lags of each variable. In VAR tests, dependencies were not consistent above the 2nd lag. Are excluded to get parsimonious computationally feasible model
                 max_lags = 2 # number lags for loop such that for all  split vars the same dataset is used (NA in the first observations otherwise)
                 lag1 = lag(base_split_variables,1)
                 lag2 = lag(base_split_variables,2)
                 colnames(lag1) = paste0(colnames(base_split_variables),"_lag1")
                 colnames(lag2) = paste0(colnames(base_split_variables),"_lag2")
                 
                 split_variables = as.data.frame(cbind(lag1 , lag2)) # data frame instead of time series

                # possible others: interactions?
                
                 
             # step 1) find optimal GARCH 1/1 for full sample. remove first max_lags obs since they are not used by TreeGarch either
                 source("scripts/garchFunction.R") # functions
                 # inputs fct
                 returns=ts_r$rub
                 ar = 1
                 ma = 1
                 threshhold = F
                 th_value  = 0 # not optimized within fct
                 data_threshhold = 0
                 type = "GARCH"
                 distribution ="normal"
                 start_parms = c(0, 0.1, 0.9,0.1)
                 opt_parms= nlm(garchEstimation,start_parms,
                                returns = returns[(max_lags+1):length(returns)],  ar = ar, ma = ma,
                                threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                                type=type, distribution=distribution,
                                print.level=2,iterlim=1000, check.analyticals=1)
               

              
          # step 2) split sample via reduction in log likelihood
               vector_quantiles = seq(1, 7)*0.125 # for threshholds
             
             # choose variable for partition
               list_split_variables = colnames(split_variables)

               split_one_covariate = drop_na(select(split_variables, list_split_variables[2])) # assign current split variable and remove first missing obs
               split_threshholds = quantile(split_one_covariate[,1],vector_quantiles)

               
             # starting values are parms of first GARCH
               start_parms = opt_parms$estimate  
              
              for (i in 1:length(split_threshholds)) {
                # split sample starting from first obs that is not NA for splitting value
                subsample1 = returns[(split_one_covariate >=split_threshholds[i])[(max_lags+1):length(returns)]] 
                subsample2 = returns[(split_one_covariate <split_threshholds[i])[(max_lags+1):length(returns)]]
                
                # estimate GARCH in subsamples and get sum of likelihood
                opt_parms1= nlm(garchEstimation,start_parms,
                               returns = subsample1,  ar = ar, ma = ma,
                               threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                               type=type, distribution=distribution,
                               iterlim=1000, check.analyticals=1)
                opt_parms2= nlm(garchEstimation,start_parms,
                                returns = subsample2,  ar = ar, ma = ma,
                                threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                                type=type, distribution=distribution,
                                iterlim=1000, check.analyticals=1)
                
                print(opt_parms$minimum)
                print(opt_parms1$minimum+ opt_parms2$minimum )
                
              }
              
            
              # step 3) prune: choose subtree that minimized sum AIC
            
           

   
          # news impact curve:
            # IMPLEMENT
               
               
          # compare results of audrino and our fct ----
          garchEstimation(theta=start_parms, returns, ar, ma, threshhold,th_value,data_threshhold,type, distribution)
          my.loglike.t( my.optpar$estimate)
          
          source("scripts/garchFunction.R") # functions
          # inputs fct
          returns=ts_r$oil
          ar = 1
          ma = 1
          threshhold = F
          th_value  = 0 # not optimized within fct
          data_threshhold = 0
          type = "GARCH"
          distribution ="normal"
          start_parms = c(0,0.1, rep(0.5,ar), rep(0.1,ma),0,6) # initialize parms
          opt_parms= nlm(garchEstimation,start_parms,
                         returns = returns,  ar = ar, ma = ma,
                         threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                         type=type, distribution=distribution,
                         print.level=1,steptol = 1e-6, iterlim=1000, check.analyticals=T)
          opt_parms$estimate
          save
          
          m1=garchFit(returns~garch(1,2),data=returns,trace=F)
          summary(m1)
          
          # Audrino fct

          par.start=c(rep(0.5,4),0,6)
          my.optpar= nlm(my.loglike.t,par.start,iterlim=1000,print.level=1)
          my.optpar$estimate
          save =  my.optpar$estimate
          
          names_coefs = c("exp_ret","constant", "MA1",  "threshhold_parm", "AR1")
          names_coefs
          c(my.optpar$estimate[1], my.optpar$estimate[2]^2,my.optpar$estimate[3]^2,my.optpar$estimate[4] ,my.optpar$estimate[5]^2 )          # make sure to revert squares in parms
          
          

          
          
          ####Univariate Garch Opt function ####
      
          # ug_spec= ugarchspec(variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(0, ), include.mean = TRUE), distribution.model="norm")
          # ug_fit= ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
          # ug_fit

          opt_garch<- function(ar,ma,returns) {
            AIC <- matrix(data=NA,nrow=2,ncol=2)
            for (j in 1:5) {
              for (i in 1:5) {
                ug_spec <- ugarchspec(variance.model=list(model="fGARCH", submodel= "TGARCH",garchOrder=c(i,j)), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="sstd") #,submodel="TGARCH"
                ugfit = ugarchfit(spec = ug_spec, data = returns, solver ='hybrid')
                info <- infocriteria(ugfit)
                AIC[j,i] <- info[2,1]
              }
            }
            opt_garch <- which(AIC == min(AIC,na.rm=TRUE), arr.ind=TRUE)
            ug_spec <- ugarchspec(variance.model=list(model="fGARCH",submodel= "TGARCH", garchOrder=c(opt_garch[1],opt_garch[2])), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="sstd") # ,submodel="TGARCH"
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
          
