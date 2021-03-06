# Data Cleaning

# libraries
  library(tidyverse)
  library(readxl)
  library(xts)
  library(quantmod)

  setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub/") # wd
  source("scripts/functions.R") # functions
  outpath =  "output/returnPlotsCleaning/"
  
# Import Data ----
  oil <- read_excel("data/Oil.xlsx")
  rub <- read_excel("data/RUBUSD.xlsx")


# Clean Data ----
  oil <- xts(oil$DCOILBRENTEU, order.by = oil$DATE)
  oil <-oil[!(oil[,1] =="."),] #NA?s are coded as "."
  
  rub <- xts(rub$`RUB/USD`, order.by = rub$data)
  
  ts <- merge(rub,oil,join='left')
  ts <- na.omit(ts)
  
  ts_r <- ts
  ts_r$rub <- dailyReturn(ts$rub, type='log')
  ts_r$oil <- dailyReturn(ts$oil, type='log')
  ts_r$rub_val <- ts$rub
  ts_r$oil_val <- ts$oil
  ts_r <- ts_r[-1, ]

# Short Descriptives: Stationarity? 
  # plot cum logreturns
    # ggplot
      title = "RUBUSD_and_Oil_Cum_Returns"
      xlab = "Time"
      ylab= "CumReturn"
      y1 = exp(cumsum(ts_r$rub))
      y2 = exp(cumsum(ts_r$oil))
      x = time(ts_r)
      names_y = c("RUB/USD", "Oil")
      line_plot_multiple(title, outpath,x,xlab, ylab, names_y, y_percent=F,y_discrete=T, legend=T, y1, y2)
      
    # normal plot
      plot(exp(cumsum(ts_r$oil)), col = 'darkgreen', main="Oil (green) and Rub/USD (blue)")
      lines(exp(cumsum(ts_r$rub)), col = 'blue',  type = 'l')
  
  # plot logreturns
      # ggplot
      title = "RUBUSD_and_Oil_Returns"
      xlab = "Time"
      ylab= "DailyReturn"
      y1 = ts_r$rub
      y2 = ts_r$oil
      x = time(ts_r)
      names_y = c("RUB/USD", "Oil")
      line_plot_multiple(title, outpath, x,xlab, ylab, names_y, y_percent=F, y_discrete=T, legend=T, y1, y2)
      
      # normal plot
      plot(ts_r$rub, main="RUB/USD")   
      plot(ts_r$oil, main="Oil")
      
  # QQ plot for Outliers and Normality (normality will be dealt with later)
      # extreme tail outliers -> remove them since no risk-mgmt application
      qqnorm(ts_r$rub)
      qqline(ts_r$rub)
      
      qqnorm(ts_r$oil)
      qqline(ts_r$oil)
      
  # Boxplots
      boxplot(coredata(ts_r$rub))
      title(main = "RUBUSD", ylab = "returns") 
      boxplot(coredata(ts_r$oil))
      title(main = "Oil", ylab = "returns") 
      
  # Remove extreme Outliers via empirical quantiles
      #quantile_outliers = 0 #disabled
      quantile_outliers = 0.01
      ts_r <- ts_r[(ts_r$rub > quantile(ts_r$rub,quantile_outliers)) & (ts_r$rub < quantile(ts_r$rub,1-quantile_outliers)) & (ts_r$oil > quantile(ts_r$oil,quantile_outliers))& (ts_r$oil < quantile(ts_r$oil,1-quantile_outliers))]
      

#save cleaned data
  #write.table(coredata(ts_r), "data/data.csv", row.names = index(ts_r), col.names = names(ts_r), sep=",")
  #write.table(coredata(ts_r), "data/data_outliers_1_with_values.csv", row.names = index(ts_r), col.names = names(ts_r), sep=",")
  ts_r <- data.frame(date=index(ts_r), coredata(ts_r))
  
  # note, that this file is overwritten by the output of datamerging.R
  write.table(ts_r, "data/data_outliers_1_with_values.csv", sep=",")
  