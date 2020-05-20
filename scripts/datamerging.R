
# libraries
library(readxl)
library(xts)
library(quantmod)
library(dplyr)

setwd("C:/Users/user/iCloudDrive/SG MiQEF/Financial Volatility/Project/GitHub/")
data <- read.table("data/data_outliers_1_with_values.csv", sep = ",")
data$date <- as.Date(data$date, "%Y-%m-%d")

## NAs are NOT dropped when merging the dataframes

# MOEX ----
data.moex <- read_excel("data/thomson_reuters/moex1.xlsx")

moex <- xts(data.moex$Close, order.by = data.moex$Date)
colnames(moex) <- c("MOEX_val")

moex$MOEX <- dailyReturn(moex$MOEX_val, type='log')
moex <- data.frame(date=index(moex), coredata(moex))
moex$date <- as.Date(moex$date, "%Y-%m-%d")
moex <- moex[-1, ]

data1 <- left_join(data, moex, by="date")

# SPX ----
data.spx <- read_excel("data/thomson_reuters/spx.xlsx")

spx <- xts(data.spx$Close, order.by = data.spx$Date)
colnames(spx) <- c("SPX_val")

spx$SPX <- dailyReturn(spx$SPX_val, type='log')
spx <- data.frame(date=index(spx), coredata(spx))
spx$date <- as.Date(spx$date, "%Y-%m-%d")
spx <- spx[-1, ]

data1 <- left_join(data1, spx, by="date")

# RUBEUR ----
data.eur <- read_excel("data/RUBEUR.xlsx")

eur <- xts(data.eur$rubeur, order.by = data.eur$date)
colnames(eur) <- c("RUBEUR_val")

eur$RUBEUR <- dailyReturn(eur$RUBEUR_val, type='log')
eur <- data.frame(date=index(eur), coredata(eur))
eur$date <- as.Date(eur$date, "%Y-%m-%d")
eur <- eur[-1, ]

data1 <- left_join(data1, eur, by="date")

# zerobond yield curve ----
data.yc <- read_excel("data/zerobond_yield_curve_ru.xlsx")

data.yc <- data.yc[, -c(3, 4)]
data.yc$date <- as.Date(data.yc$date, "%Y-%m-%d")

data1 <- left_join(data1, data.yc, by="date")

# 3m tbill ----
data.tt <- read.table("data/fred/3m_tbill.csv", sep=",", header = TRUE)
data.tt <-data.tt[!(data.tt[,2] =="."),]

data.tt$DATE <- as.Date(data.tt$DATE, "%Y-%m-%d")
colnames(data.tt) <- c("date", "DGS3MO")

data1 <- left_join(data1, data.tt, by="date")

# 10 tbill ----
data.tet <- read.table("data/fred/10y_tbill.csv", sep=",", header = TRUE)
data.tet <-data.tet[!(data.tet[,2] =="."),]

data.tet$DATE <- as.Date(data.tet$DATE, "%Y-%m-%d")
colnames(data.tet) <- c("date", "DGS10")

data1 <- left_join(data1, data.tet, by="date")

# save the data
write.table(data1, "data/data_outliers_1_with_values.csv", sep=",")
