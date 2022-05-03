# Dependence of Exchange Rate Volatility on Commodity Prices Example of Russian Ruble

# Project
This project examines the dependence of exchange rate volatility on commodity prices.
We analyse the example of the Russian Ruble / US-Dollar exchange rate and its dependence on crude oil prices.
A simple example of a tree-structured GARCH model for volatility is provided including our own implementation in R.

# Paper Abstract
This paper examines the relationship between volatility of of the Russian
Ruble / US Dollar exchange rate and volatility of the Brent Oil
price. For this we use a classic GARCH model, a tree-GARCH model and
a DCC-GARCH model. We also test for structural breaks. We show that
at times of decreasing oil price returns, the volatility of the Russian Ruble
/ US Dollar exchange rate is less persistent and the conditional correlation
between Russian Ruble / US Dollar exchange rate returns and Brent oil
price returns is more persistent than in normal market conditions.

# Contributors
Johannes Cordier, Liudmila Gorkun-Voevoda, Erik-Jan Senn

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
