# Volatile-Gamma

# Data sources
Exchange rate:
https://www.cbr.ru/currency_base/dynamics/?UniDbQuery.Posted=True&UniDbQuery.mode=1&UniDbQuery.date_req1=&UniDbQuery.date_req2=&UniDbQuery.VAL_NM_RQ=R01235&UniDbQuery.From=01.01.2000&UniDbQuery.To=15.05.2020

Oil:
https://fred.stlouisfed.org/series/DCOILBRENTEU

# Work distribution
* Erik
  * Test tree GARCH (univariate)
  * Write a GARCH function
  * Granger causality
  * look for more data
  * research question

* Mila
  * Cleaning file
  * AutoARIMA -> save parameters, give the residuals to GARCH estimation (how to chosse tgarch threshold?)
  * Look at tail distributions (t, GED only if possible in multivariate garch
  
* Johannes
  * DCC
  * Try to keep flexible the number of univariate GARCHes
