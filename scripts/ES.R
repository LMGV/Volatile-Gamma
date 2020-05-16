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
  library(rugarch)
  library(rmgarch)
  library(psych)
  library(MASS)
  library(fitdistrplus)
  
  filter <- dplyr::filter
  select <- dplyr::select

source("scripts/functions.R") # functions

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
          
          ####Univariate Garch Opt function ####
          opt_garch<- function(ar,ma,returns) {
            AIC <- matrix(data=NA,nrow=2,ncol=2)
            for (j in 1:2) {
              for (i in 1:2) {
                ug_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(j,i)), mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE), distribution.model="sstd") #,submodel="TGARCH"
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
