# Data Cleaning

# librarys
  library(readxl)
  library(xts)
  library(quantmod)

  setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub/data")

# Import Data ----
  oil <- read_excel("Oil.xlsx")
  rub <- read_excel("RUBUSD.xlsx")

# Clean Data ----
  oil <- xts(oil$DCOILBRENTEU, order.by = oil$DATE)
  oil<-oil[!(oil[,1] =="."),] #NA?s are coded as "."
  
  rub <- xts(rub$`RUB/USD`, order.by = rub$data)
  
  ts <- merge(rub,oil,join='left')
  ts <- na.omit(ts)
  
  ts_r <- ts
  ts_r$rub <- dailyReturn(ts$rub, type='log')
  ts_r$oil <- dailyReturn(ts$oil, type='log')
  ts_r <- ts_r[-1, ]

# Short Descriptives: Stationarity? 
  # plot cum logreturns
  plot(exp(cumsum(ts_r$oil)), col = 'darkgreen', main="Oil (green) and Rub/USD (blue)")
  lines(exp(cumsum(ts_r$rub)), col = 'blue',  type = 'l')
  
  # plot logreturns
  plot(ts_r$rub)
  plot(ts_r$oil)

#save cleaned data
  write.table(coredata(ts_r), "data.csv", row.names = index(ts_r), col.names = names(ts_r), sep=",")