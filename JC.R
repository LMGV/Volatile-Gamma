library("readxl")    
library("ggplot2")    
library("ggfortify")  
library("zoo")        
library("xtable")     
library("lmtest")     
library("forecast") 
library(tseries)
library(xts)
library(quantmod)

setwd("C:/Users/johan/Google Drive/Vola")


####Data Cleaning####
oil <- read_excel("Oil.xlsx")
rub <- read_excel("RUBUSD.xlsx")

oil <- xts(oil$DCOILBRENTEU, order.by = oil$DATE)
oil<-oil[!(oil[,1] =="."),] #NA´s are coded as "."

rub <- xts(rub$`RUB/USD`, order.by = rub$data)

ts <- merge(rub,oil,join='left')
ts <- na.omit(ts)

ts_r <- ts
ts_r$rub <- dailyReturn(ts$rub, type='log')
ts_r$oil <- dailyReturn(ts$oil, type='log')


####Stationarity####
plot(ts_r$rub)
plot(ts_r$oil)

####Mean####
print("Mean of Rubel log returns")
mean(ts_r$rub)

print("Mean of oil log returns")
mean(ts_r$oil)


####ACF and PACF of r####
acf(ts_r$oil, lag.max = 30, plot = TRUE)
acf(ts_r$rub, lag.max = 30, plot = TRUE)


####ACF and PACF of r^2####
ts_r2 <- ts_r^2
acf(ts_r2$oil, lag.max = 30, plot = TRUE)
acf(ts_r2$rub, lag.max = 30, plot = TRUE)

####ARIMA####
auto.arima(ts_r$oil)
auto.arima(ts_r$rub)
