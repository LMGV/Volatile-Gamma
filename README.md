# Volatile-Gamma

# Overleaf

https://www.overleaf.com/8847323888ghkkzmmqtgwb 


# Data sources
Exchange rate:
https://www.cbr.ru/currency_base/dynamics/?UniDbQuery.Posted=True&UniDbQuery.mode=1&UniDbQuery.date_req1=&UniDbQuery.date_req2=&UniDbQuery.VAL_NM_RQ=R01235&UniDbQuery.From=01.01.2000&UniDbQuery.To=15.05.2020

Oil & riskfree rates:
https://fred.stlouisfed.org/series/DCOILBRENTEU
https://fred.stlouisfed.org/series/IRLTCT01RUM156N
https://fred.stlouisfed.org/series/DGS10
https://fred.stlouisfed.org/series/DGS3MO

Riskfree rates Russian Federation
https://cbr.ru/hd_base/zcyc_params/

Russian Government Bond Zero Coupon Yield Curve:
https://cbr.ru/eng/hd_base/zcyc_params/?UniDbQuery.Posted=True&UniDbQuery.From=14.05.2010&UniDbQuery.To=15.05.2020

Thomson Reuters for Exchange rates and Indices 

# Work distribution
* Erik
  * Test tree GARCH (univariate)
  * Write a GARCH function
  * Granger causality
  * look for more data
	SPX, MOEX (russian index, denoted in RUB, 47% energy sector contribution), US 3m and 10y rf.
	Missing: RU 3m and 10y rf
  * research question
  * what to do with weekends / other missing values in one series? 
	Result: common consensus: ignore weekends. drop days where not all observations for all times series are given (on return level)

* Mila
  * Cleaning file
  * AutoARIMA -> save parameters, give the residuals to GARCH estimation (how to chosse tgarch threshold?)
  * Look at tail distributions (t, GED only if possible in multivariate garch
  
* Johannes
  * DCC
  * Try to keep flexible the number of univariate GARCHes
