# Data Cleaning

# librarys
  library(tidyverse)
  library(readxl)
  library(xts)
  library(quantmod)

  setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub/") # wd
  source("scripts/functions.R") # functions
  
# Import Data ----
  oil <- read_excel("data/Oil.xlsx")
  rub <- read_excel("data/RUBUSD.xlsx")


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
    # ggplot
      title = "Oil_and_RUBUSD_Cum_Returns"
      xlab = "Time"
      ylab= "CumReturn"
      y1 = exp(cumsum(ts_r$oil))
      y2 = exp(cumsum(ts_r$rub))
      x = time(ts_r)
      names_y = c("Oil","RUB/USD")
      line_point_plot_multiple(title, x,xlab, ylab, names_y, y_percent=F,y_discrete=T, legend=T, y1, y2)
      
    # normal plot
      plot(exp(cumsum(ts_r$oil)), col = 'darkgreen', main="Oil (green) and Rub/USD (blue)")
      lines(exp(cumsum(ts_r$rub)), col = 'blue',  type = 'l')
  
  # plot logreturns
      # ggplot
      title = "Oil_and_RUBUSD_Returns"
      xlab = "Time"
      ylab= "DailyReturn"
      y1 = ts_r$rub
      y2 = ts_r$oil
      x = time(ts_r)
      names_y = c("Oil","RUB/USD")
      line_point_plot_multiple(title, x,xlab, ylab, names_y, y_percent=F, y_discrete=T, legend=T, y1, y2)
      
      # normal plot
      plot(ts_r$oil, main="Oil")
      plot(ts_r$rub, main="RUB/USD")


#save cleaned data
  write.table(coredata(ts_r), "data/data.csv", row.names = index(ts_r), col.names = names(ts_r), sep=",")