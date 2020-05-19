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
  * 
  * write garchPredict fct
  * 
  * GARCH model analysis: 
	ACF / PACF + PF tests, updated ConfIntervalls for ARIMIA. Investigate before and after model selection

* Mila
  * Structural Breaks: see file arimaStructural Breaks Start with comments @Mila where you need to change things
	idea: LR test for struc break, grid of breaks each quarter
	identify best struc break and then do test again
	predict based on selected arimas
	save this again as ONE file for entire timeframe
  * AutoARIMA -> save parameters, give the residuals to GARCH estimation (how to chosse tgarch threshold?)
  * Look at tail distributions (t, GED only if possible in multivariate garch
  
* Johannes
  * DCC
  * Try to keep flexible the number of univariate GARCHes
  
# Ideas next meeting
* Structure / Code: 	
	flexible paths (setwd("~/GitHub/Volatile-Gamma")), for inputs and outputs
	always save same type of file as in cleaning
	
* likelihood line: GARCH 1/1 normal : line 506, 634(t), 755

* Test for structural breaks (different coefs within time). On uni or multivariate level?. LR Test
* Mincer-Zarnowitz-Regressions for univariate RUB. maybe for dcc ->evaluated model performance

# Work done / Cancelled
* Erik
  * Write a GARCH function
  * Granger causality - postponed: Vector Autoregressive Process, need to test for sign parms. Need to check specifics
  * look for more data - done
	SPX, MOEX (russian index, denoted in RUB, 47% energy sector contribution), US 3m and 10y rf.
	Missing: RU 3m and 10y rf
  * research question
  * HELP: likelihoods for: GARCH (2,1), GARCH(1,2), GARCH (2,2) for Normal and T-distrib. Understand likelihood in audrino-code)
  * generalize GARCH fct for 2/3 lags?
	what to do with weekends / other missing values in one series? - done (just remove, its common practice)
	Result: common consensus: ignore weekends. drop days where not all observations for all times series are given (on return level)
* Mila
  * Cleaning file
* Johannes
