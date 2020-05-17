##R-Code for the lecture "Financial Volatility"
################################################

(most data used in the examples can be found on the website

http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts3/

)

setwd("D:/Documenti/Lezioni/Financial Volatility/Notes")

#Chapter 1, page 9:

library(fBasics) #Load the package fBasics
#ibm.sp500=read.table("d-ibm3dx7008.txt",header=T) #Load the data IBM, S&amp;P500
#intel=read.table("d-intc7208.txt",header=T) #Load the data Intel
#microsoft=read.table("d-msft8608.txt",header=T) #Load the data Microsoft

basicStats(ibm.sp500[,2]*100) #Compute the summary statistics
#Alternative: use commands mean(), var(), skewness(), kurtosis()

#Simple tests:
s1=skewness(ibm.sp500[,2]*100)
t1=s1/sqrt(6/9845)
2*(1-pnorm(t1)) #Compute p-value

#Turn to log returns in percentages:
libm=log(ibm.sp500[,2]+1)*100
t.test(libm) #Test mean being zero

normalTest(libm,method="jb") #Jarque-Bera test for normality

#Chapter 1, page 12:
a<-seq(-4,4,0.01)
plot(a,dnorm(a),type="l",xlab="x",ylab="f(x)")
lines(a,dcauchy(a),lty=2)
lines(a,(1-0.10)*dnorm(a)+0.10*dnorm(a,sd=4),lty=3)

#Chapter 1, page 20:

#ibm.sp500.m=read.table("m-ibm3dx2608.txt",header=T) #Load the data IBM, S&amp;P500 monthly
libm=log(ibm.p500.m[,2]+1)
lsp500=log(ibm.sp500.m[,4]+1)
par(mfrow=c(2,1))
acf(libm,lag.max=50,main="IBM")
acf(lsp500,lag.max=50,main="S&P500")

Box.test(lsp500,lag=5,type="Ljung") #Ljung-Box statistic Q(5)

#Chapter 1, pages 24/25:

par(mfrow=c(2,2))
set.seed(66)      # so you can reproduce these results
a = arima.sim(list(order=c(1,0,0), ar=.9), n=100) #AR(1)  
acf(a,lag.max=30,main="AR(1)")
a = arima.sim(list(order=c(1,0,0), ar=-0.8), n=100) #AR(1)  
acf(a,lag.max=30,main="AR(1)")
a = arima.sim(list(order=c(2,0,0), ar=c(1.2,-0.35)), n=100) #AR(2)
acf(a,lag.max=30,main="AR(2)")
a = arima.sim(list(order=c(2,0,0), ar=c(-0.2,0.35)), n=100) #AR(2)
acf(a,lag.max=30,main="AR(2)")

par(mfrow=c(2,2))
set.seed(148)
a = arima.sim(list(order=c(0,0,1), ma=-0.7), n=100) #MA(1)
acf(a,lag.max=30,main="MA(1)")
a = arima.sim(list(order=c(0,0,2), ma=c(-0.7,0.7)), n=100) #MA(2)
acf(a,lag.max=30,main="MA(2)")
a = arima.sim(list(order=c(1,0,1), ma=0.1,ar=0.9), n=100) #ARMA(1,1)
acf(a,lag.max=30,main="ARMA(1,1)")
a = arima.sim(list(order=c(1,0,1), ma=-0.77,ar=0.55), n=100) #ARMA(1,1)
acf(a,lag.max=30,main="ARMA(1,1)")

#For autoregressive models use PACF to have an idea of the order:
a = arima.sim(list(order=c(1,0,0), ar=-0.8), n=100) #AR(1)
a = arima.sim(list(order=c(2,0,0), ar=c(1.2,-0.35)), n=100) #AR(2)

pacf(a)

#Chapter 1, page 27:

#ew.crsp=read.table("m-ew6299.txt",header=F) #Load the data

par(mfrow=c(1,2))
acf(ew.crsp,main="CRSP equal-weighted index")
pacf(ew.crsp,main="CRSP equal-weighted index")

ew.crsp=ts(ew.crsp,frequency=12,start=c(1962,1)) #Time series transformation
ts.plot(ew.crsp)

a=ar(ew.crsp,max.order=10) #Optimal order
a$order
a=arima(ew.crsp,order=c(1,0,0)) #Estimation
a=arima(ew.crsp,order=c(0,0,1)) #Optimal order MA(1) by graphical inspection of ACF


#Chapter 1, pages 31-33:

#gdp=read.table("q-gdpdef.txt",header=T) #Load the data
#gdp=ts(gdp[,4],frequency=4,start=c(1947,1))

#Step (i):
par(mfrow=c(3,1))
ts.plot(gdp,main="GDP deflator",ylab="")

gdp.fd=ts(gdp[2:249]-gdp[1:248],frequency=4,start=c(1947,2))
ts.plot(gdp.fd,main="GDP deflator first differences",ylab="") #First differences

gdp.sd=ts(gdp[3:249]-2*gdp[2:248]+gdp[1:247],frequency=4,start=c(1947,3))
ts.plot(gdp.sd,main="GDP deflator second differences",ylab="") #Second differences

#Step (ii):
par(mfrow=c(1,2))
acf(gdp.fd[-248],main="GDP deflator first differences")
pacf(gdp.fd[-248],main="GDP deflator first differences")

#Step (iii):
a=arima(gdp.fd[-248],order=c(4,0,0))

#Step (iv):
par(mfrow=c(1,2))
acf(a$residuals,main="AR(4) residuals")
pacf(a$residuals,main="AR(4) residuals")

Box.test(a$residuals,lag=50,type="Ljung") #Ljung-Box statistic Q(5)Box.test

normalTest(a$residuals,method="jb") #Jarque-Bera test

#Chapter 1, pages 35ss:

library(zoo)
a=strptime(ibm.sp500[,1],"%Y%m%d")
libm=log(ibm.sp500[,2]+1)*100
lsp500=log(ibm.sp500[,4]+1)*100

par(mfrow=c(2,1))
plot(zoo(libm,a),ylab="IBM returns",xlab="")
plot(zoo(lsp500,a),ylab="S&P500 returns",xlab="")

par(mfrow=c(1,2))
acf(libm,lag.max=50,main="IBM: returns")
acf(libm^2,lag.max=50,main="IBM: squared returns")

a=c()
for (i in 1:length(libm))
{
  a=c(a,max(libm[i],0))
}

for (h in 1:7)
{
  print(cor(a[1:(length(a)-h)],abs(libm)[(1+h):length(libm)])) #Leverage effect
}

a=c()
for (i in 1:length(libm))
{
  a=c(a,max(-libm[i],0))
}

for (h in 1:7)
{
  print(cor(a[1:(length(a)-h)],abs(libm)[(1+h):length(libm)])) #Leverage effect
}


par(mfrow=c(1,2))
libm=log(ibm.sp500[,2]+1)*100
a=density(libm)
plot(a,main="Daily IBM returns kernel density",xlim=c(-10,10))
lines(a$x,dnorm(a$x,mean=mean(libm),sd=sqrt(var(libm))),lty=2)

libm=log(ibm.sp500.m[,2]+1)*100
a=density(libm)
plot(a,main="Monthly IBM returns kernel density",xlim=c(-40,40))
lines(a$x,dnorm(a$x,mean=mean(libm),sd=sqrt(var(libm))),lty=2)

##########################################################################

#Chapter 2, page 8ss:

intel.m=read.table("m-intc7308.txt",header=T) #Load the data Intel
intel.m=log(intel.m[,2]+1)
intel.m=ts(intel.m,frequency=12,start=c(1973,1))

par(mfrow=c(2,2))
plot(intel.m,main="Intel monthly log returns")
acf(intel.m,main="Intel log returns")
acf(intel.m^2,main="Intel squared returns")
pacf(intel.m^2,main="Intel squared returns")

Box.test(intel.m,lag=30,type="Ljung")
Box.test(intel.m^2,lag=30,type="Ljung")

library(fGarch) #Load the package for ARCH/GARCH estimation
m1=garchFit(intel.m~garch(3,0),data=intel.m,trace=F)
summary(m1) #Obtain results

m1=garchFit(intel.m~garch(1,0),data=intel.m,trace=F)
summary(m1) #Obtain results
predict(m1,5) #Predictions until month t+5

m1=garchFit(intel.m~garch(1,0),data=intel.m,trace=F,cond.dist="std")
summary(m1)

par(mfrow=c(2,2))
qqplot(qt(ppoints(length(m1@residuals)), df = 6),y=m1@residuals/m1@sigma.t ,main = expression("Q-Q plot for" ~~ {t}[nu == 6]),xlab="qt",ylab="residuals")
acf(m1@residuals/m1@sigma.t,main="ARCH(1) residuals")
acf(m1@residuals^2/m1@sigma.t^2,main="ARCH(1) squared residuals")
acf(abs(m1@residuals/m1@sigma.t),main="ARCH(1) absolute residuals")

#Extensions:
#m1=garchFit(intel.m~arma(1,0)+garch(1,0),data=intel.m,trace=F)

#Chapter 2, page 14ss:

par(mfrow=c(2,1))
set.seed(122)
spec = garchSpec(model = list(omega= 1, alpha = 0.95, beta = 0)) #ARCH(1) model
a=garchSim(spec, n = 500, n.start=100)
ts.plot(unclass(a),main="",xlab="",ylab="")

par(mfrow=c(2,1))
set.seed(333)
spec = garchSpec(model = list(omega= 1, alpha = 0.7, beta = 0.2)) #GARCH(1,1) model
a=garchSim(spec, n = 500, n.start=100)
ts.plot(unclass(a),main="",xlab="",ylab="")

#Chapter 2, page 25:

par(mfrow=c(3,2))
ibm=log(ibm.sp500[,2]+1)
microsoft=log(microsoft.dat[,2]+1)
sp500=log(ibm.sp500[,4]+1)
acf(ibm^2,main="IBM")
pacf(ibm^2,main="IBM")
acf(microsoft^2,main="Microsoft")
pacf(microsoft^2,main="Microsoft")
acf(sp500^2,main="S&amp;P500")
pacf(sp500^2,main="S&amp;P500")

#Chapter 2, page 27:

par(mfrow=c(2,1))

spec = garchSpec(model = list(omega= 1, alpha = 0.1, beta = 0.8)) #GARCH(1,1) model

pred.for<-c()
for (t in 1:500)
{
   set.seed(56)
   a=garchSim(spec, n = 500+t-1, n.start=100)
   a=garchFit(formula = ~ garch(1, 1), data = a,include.mean=FALSE,trace=FALSE)
   b=predict(a,n.ahead=1)
   print(b)
   pred.for<-c(pred.for, qnorm(0.975)*b$standardDeviation)
}

set.seed(56)
a=garchSim(spec, n = 1000, n.start=100)
ts.plot(a$garch[501:1000],ylim=c(-11,11),xlab="",ylab="")
lines(pred.for,lty=2)
lines(-pred.for,lty=2)

spec = garchSpec(model = list(omega= 1, alpha = 0.6, beta = 0.2)) #GARCH(1,1) model

pred.for<-c()
for (t in 1:500)
{
   set.seed(224)
   a=garchSim(spec, n = 500+t-1, n.start=100)
   a=garchFit(formula = ~ garch(1, 1), data = a,include.mean=FALSE,trace=FALSE)
   b=predict(a,n.ahead=1)
   print(b)
   pred.for<-c(pred.for, qnorm(0.975)*b$standardDeviation)
}

set.seed(224)
a=garchSim(spec, n = 1000, n.start=100)
ts.plot(a$garch[501:1000],ylim=c(-20,20),xlab="",ylab="")
lines(pred.for,lty=2)
lines(-pred.for,lty=2)

#Chapter 2, page 28:

m1=garchFit(intel.m~garch(1,1),data=intel.m,trace=F)
summary(m1) #Obtain results

par(mfrow=c(1,2))
qqnorm(m1@residuals/m1@sigma.t)
#acf(m1@residuals/m1@sigma.t,main="GARCH(1,1) residuals")
acf(m1@residuals^2/m1@sigma.t^2,main="GARCH(1,1) squared residuals")
#acf(abs(m1@residuals/m1@sigma.t),main="GARCH(1,1) absolute residuals")



#Chapter 2, page 35:

gamma=function(x,h)
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  gamma=sum(x[1:(n-h)]*x[(h+1):n])/n
}

rho=function(x,h)
{
  rho=gamma(x,h)/gamma(x,0)
}

n1.acf=function(x,main=NULL,method="NP")
{
  n=length(x)
  nlag=as.integer(min(10*log10(n),n-1))
  acf.val=sapply(c(1:nlag),function(h) rho(x,h))
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(acf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(acf.val,1.96*band,1.96/sqrt(n))
  acf(x,xlab="Lag",ylab="Sample autocorrelations",ylim=c(minval,maxval),main=main)
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}

par(mfrow=c(2,2)) 
n1.acf(microsoft,main=c("Microsoft")) #daily data
n1.acf(ibm,main=c("IBM")) #daily data
n1.acf(unclass(intel.m)[1:length(intel.m)],main=c("Intel")) #monthly data
n1.acf(sp500,main=c("S&amp;P500")) #daily data

#Chapter 2, page 36:

Box.test(microsoft,lag=20,type="Ljung") #Ljung-Box statistic Q(5) Box.test

gamma.asy=function(x,h) #Asymptotical cov. matrix of the sample autocorrelations
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  x2=x^2
  gamma.asy<-matrix(,h,h)  
  for (i in 1:h)
  {
    for (j in i:h)
    {
      gamma.asy[i,j]=gamma.asy[j,i]=sum(x[(j-i+1):(n-i)]*x[1:(n-j)]*x2[(j+1):n])/n
    }
  }
  rho.asy=1/gamma(x,0)^2*gamma.asy
  list(gamma.asy=gamma.asy,rho.asy=rho.asy)
}

#a=gamma.asy(microsoft,5)
corr.Box.test=function(x,h) #Corrected portmanteau test under GARCH assumption
{
  n<-length(x)
  a=gamma.asy(x,h)
  acf.val=sapply(c(1:h),function(h) rho(x,h))
  val=n*(acf.val%*%solve(a$rho.asy)%*%acf.val)
  print(val)  
  print(1-pchisq(val,h))
}

corr.Box.test(microsoft,20)

#Chapter 2, page 42/43:
Box.test(sp500^2,lag=10,type="Ljung")

LM=function(x,h)
{
  n=length(x)
  x2=x^2-mean(x^2)
  dat=matrix(,n-h,h+1)
  for (i in 1:(h+1))
   {
     dat[,i]=x2[(h+2-i):(n-i+1)]
   }
  a=lm(dat[,1]~dat[,2:(h+1)])
  r2=summary(a)$r.squared
  print(r2 * n)
  print(1-pchisq(r2*n,h))
}

LM(sp500,10)

par(mfrow=c(2,2))
acf(microsoft.m^2,main="Microsoft")
acf(ibm.m^2,main="IBM")
acf(intel.m^2,main="Intel")
acf(sp500.m^2,main="sp500")

#Chapter 2, page 50:

library(fGarch) #Load the package for ARCH/GARCH estimation
m1=garchFit(microsoft~garch(1,1),data=microsoft,trace=F)
summary(m1) #Obtain results

(m1@fit$matcoef[3,1]+m1@fit$matcoef[4,1])
(m1@fit$matcoef[3,1]+m1@fit$matcoef[4,1])^2+(mean(m1@residuals^4/m1@sigma.t^4)-1)*m1@fit$matcoef[3,1]^2  #rho4

#Chapter 2, page 52ss:
m1=garchFit(microsoft~garch(1,1),data=microsoft,trace=F)
summary(m1) #Obtain results

par(mfrow=c(2,1))
acf(microsoft^2)
pacf(microsoft^2)

library(zoo)

#microsoft.dat=read.table("d-msft8608.txt",header=T)
a=strptime(microsoft.dat[,1],"%Y%m%d")

par(mfrow=c(3,1))
plot(zoo(microsoft,a),ylab="",main="Microsoft daily returns",xlab="")

m1=garchFit(microsoft~garch(1,1),data=microsoft,trace=F)
plot(zoo(m1@sigma.t,a),xlab="",ylab="",main="Conditional volatility from GARCH(1,1)")

m1=garchFit(microsoft~garch(5,0),data=microsoft,trace=F)
plot(zoo(m1@sigma.t,a),xlab="",ylab="",main="Conditional volatility from ARCH(5)")

par(mfrow=c(2,2))
plot(zoo(m1@residuals/m1@sigma.t,a),xlab="",ylab="",main="GARCH(1,1) residuals")
acf(m1@residuals/m1@sigma.t,main="GARCH(1,1) residuals")
acf(m1@residuals^2/m1@sigma.t^2,main="GARCH(1,1) squared residuals")
qqnorm(m1@residuals/m1@sigma.t)

#######################################################################

#Chapter 3, page 5:

a=c()
for (i in 1:length(ibm))
{
  a=c(a,max(sp500[i],0))
}

for (h in 1:40)
{
  print(h)
  print(cor(a[(1+h):(length(a))],ibm[(1):(length(sp500)-h)])) #Leverage effect
}

a=(sp500<0) #Sign Bias
a=c()
for (i in 1:length(ibm))
{
  a=c(a,100*min(sp500[i],0)) #Negative Size Bias
  #a=c(a,100*max(sp500[i],0)) #Positive Size Bias
}

h=5
b=lm(100*sp500[(h+1):length(ibm)]^2~a[1:(length(ibm)-h)])
summary(b)

#Chapter 3, page 18:

a=seq(-10,10,0.01)
plot(a,1+0.38*a^2,type="l",lty=1,main="News impact curve",xlab="Past innovation",ylab="Conditional variance") #ARCH(1)
lines(a,(1-0.5*a*(a<=0)+0.2*a*(a>0))^2,lty=2,col="blue") #TARCH(1)
lines(a,1+0.2*a^2+0.36*a^2*(a<=0),lty=3,col="darkred") #GJR-ARCH(1)
lines(a,1+0.2*(abs(a)-0.948*a)^2,lty=4,col="darkgreen") #PARCH(1) with d=2

#Chapter 3, page 19ss:

my.loglike.normal=function(theta) #Estimate an asymmetric GARCH(1,1) model
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,n))
  log.sigmasq=c(log(my.sigmasq[1]),rep(0,n))
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(n+1))
  {
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    #aiuto=(data[i-1]-my.mean[i-1])/sqrt(exp(log.sigmasq[i-1]))
    #log.sigmasq[i]=theta[2]+theta[3]*(theta[4]*aiuto+abs(aiuto)-sqrt(2/pi))+theta[5]*log.sigmasq[i-1] #EGARCH(1,1) 
    
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)
  1/2*sum(log(my.sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]))))
 
}

returns=ibm
par.start=rep(0.5,5)
my.optpar= nlm(my.loglike.normal,par.start,iterlim=1000,print.level=2)


sigmasq.model=function(theta)
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=sqrt(my.sigmasq)
  log.sigmasq=log(my.sigmasq)
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(n+1))
  {
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    #aiuto=(data[i-1]-my.mean[i-1])/sqrt(exp(log.sigmasq[i-1]))
    #log.sigmasq[i]=theta[2]+theta[3]*(theta[4]*aiuto+abs(aiuto)-sqrt(2/pi))+theta[5]*log.sigmasq[i-1] #EGARCH(1,1) 
    
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)


	list(my.sigmasq = my.sigmasq[2:(n + 1)],my.mean=my.mean[2:(n+1)])
}

a=sigmasq.model(my.optpar$estimate)

#Save estimates and optimal parameters in ibm.garch11, ibm.egarch11, ...

#Estimates news impact curves

aiuto=var(ibm)
a=seq(min(ibm-mean(ibm))-0.1,max(ibm-mean(ibm))+0.2,0.001)
#aiuto=2.292193e-06/(1-6.160892e-02-9.332932e-01)
plot(a,2.292193e-06+6.160892e-02*a^2+9.332932e-01*aiuto,type="l",lty=1,main="Estimated news impact curves for the IBM data",xlab="Past innovation",ylab="Conditional variance",ylim=c(0,0.01)) #GARCH(1,1)
#aiuto=2.308737e-06/(1-2.530091e-02-6.851351e-02/2-9.360238e-01)
lines(a,2.308737e-06+2.530091e-02*a^2+6.851351e-02*a^2*(a<=0)+9.360238e-01*aiuto,lty=3,col="blue") #GJR-GARCH(1,1)
#aiuto=0.0001149456^2/(1-0.0610985383/sqrt(2/pi)-0.9472648804)^2
A=(0.0001149456+0.9472648804*sqrt(aiuto))^2 	
lines(a,A+2*sqrt(A)*0.0610985383*(abs(a)-0.4603728320*a)+0.0610985383^2*(abs(a)-0.4603728320*a)^2,lty=2,col="darkred") #TGARCH(1,1)


library(zoo)
a=strptime(ibm.sp500[,1],"%Y%m%d")

par(mfrow=c(2,2))
plot(zoo(ibm.garch11$my.sigmasq,a),main="GARCH",xlab="",ylab="")
plot(zoo(ibm.egarch11$my.sigmasq,a),main="EGARCH",xlab="",ylab="")
plot(zoo(ibm.gjr11$my.sigmasq,a),main="GJR-GARCH",xlab="",ylab="")
plot(zoo(ibm.tgarch11$my.sigmasq,a),main="TGARCH",xlab="",ylab="")

par(mfrow=c(1,2))
plot(ibm,main="GARCH",xlab="",ylab="",type="l",ylim=c(-0.25,0.25))
lines(mean(ibm)+3*sqrt(ibm.garch11$my.sigmasq),lty=4,col="blue")
lines(mean(ibm)-3*sqrt(ibm.garch11$my.sigmasq),lty=4,col="blue")

plot(ibm,main="TGARCH",xlab="",ylab="",type="l",ylim=c(-0.25,0.25))
lines(mean(ibm)+3*sqrt(ibm.tgarch11$my.sigmasq),lty=4,col="blue")
lines(mean(ibm)-3*sqrt(ibm.tgarch11$my.sigmasq),lty=4,col="blue")

sum((ibm>(mean(ibm)+3*sqrt(ibm.gjr11$my.sigmasq))))+sum((ibm<(mean(ibm)-3*sqrt(ibm.gjr11$my.sigmasq)))) #Compute mumber of violations

mean((sqrt(ibm.tgarch11$my.sigmasq)-sqrt(ibm.egarch11$my.sigmasq))^2) #Distance between volatility estimates

plot(ibm.tgarch11$my.sigmasq,ibm.egarch11$my.sigmasq,main="TGARCH vs. EGARCH",xlab="",ylab="")
plot(ibm.tgarch11$my.sigmasq,ibm.gjr11$my.sigmasq,main="TGARCH vs. GJR-GARCH",xlab="",ylab="")

#Chapter 3, page 26ss:

my.loglike.t=function(theta) #Estimate an asymmetric GARCH(1,1) model with Student's t innovations
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,n))
   
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  for (i in 2:(n+1))
  {
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)
  
  1/2*sum(log(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
 }

returns=ibm
par.start=c(rep(0.5,5),4)
my.optpar= nlm(my.loglike.t,par.start,iterlim=1000,print.level=2)

sigmasq.model=function(theta)
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=sqrt(my.sigmasq)
 
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(n+1))
  {
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)


	list(my.sigmasq = my.sigmasq[2:(n + 1)],my.mean=my.mean[2:(n+1)])
}

a=sigmasq.model(my.optpar$estimate)


library(zoo)
a=strptime(ibm.sp500[,1],"%Y%m%d")

par(mfrow=c(2,3))
plot(zoo(ibm.garch11$my.sigmasq,a),main="GARCH",xlab="",ylab="")
plot(zoo(ibm.gjr11$my.sigmasq,a),main="GJR-GARCH",xlab="",ylab="")
plot(zoo(ibm.tgarch11$my.sigmasq,a),main="TGARCH",xlab="",ylab="")
plot(zoo(ibm.garch11.t$my.sigmasq,a),main="GARCH-t",xlab="",ylab="")
plot(zoo(ibm.gjr11.t$my.sigmasq,a),main="GJR-GARCH-t",xlab="",ylab="")
plot(zoo(ibm.tgarch11.t$my.sigmasq,a),main="TGARCH-t",xlab="",ylab="")

mean((sqrt(ibm.gjr11.t$my.sigmasq)-sqrt(ibm.garch11$my.sigmasq))^2) #Distance between volatility estimates

par(mfrow=c(1,2))
plot(ibm.garch11$my.sigmasq,ibm.garch11.t$my.sigmasq,main="GARCH vs. GARCH-t",xlab="",ylab="")
plot(ibm.tgarch11.t$my.sigmasq,ibm.gjr11.t$my.sigmasq,main="TGARCH-t vs. GJR-GARCH-t",xlab="",ylab="")

par(mfrow=c(2,2))
qqnorm((returns-ibm.garch11$my.mean)/sqrt(ibm.garch11$my.sigmasq),main="GARCH")
qqplot((returns-ibm.garch11.t$my.mean)/sqrt(ibm.garch11.t$my.sigmasq), rt(length(returns), df = ibm.garch11.t$opt.par[5]),main="GARCH-t",xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqnorm((returns-ibm.tgarch11$my.mean)/sqrt(ibm.tgarch11$my.sigmasq),main="TGARCH")
qqplot((returns-ibm.tgarch11.t$my.mean)/sqrt(ibm.tgarch11.t$my.sigmasq),rt(length(returns), df = ibm.tgarch11.t$opt.par[6]),main="TGARCH-t",xlab="Theoretical Quantiles",ylab="Sample Quantiles")

#library(tseries)
#jarque.bera.test((ibm-ibm.garch11$my.mean)/sqrt(ibm.garch11$my.sigmasq))

sum((ibm>(mean(ibm)-qt(0.001,df=6.2)*sqrt(ibm.gjr11.t$my.sigmasq))))+sum((ibm<(mean(ibm)+qt(0.001,df=6.2)*sqrt(ibm.gjr11$my.sigmasq)))) #Compute mumber of violations

#Chapter 3, page 30:
par(mfrow=c(2,2))
acf(ibm^2,lag.max=50,main="IBM")
acf(intel^2,lag.max=50,main="Intel")
acf(microsoft^2,lag.max=50,main="Microsoft")
acf(sp500^2,lag.max=50,main="S&P500")

#Chapter 3, page 33:

a=c()
for (i in 1:length(ibm))
{
  a=c(a,sum(ibm[1:i]^2-mean(ibm^2)))
}

1/sqrt(length(ibm))*1/sqrt(var(ibm^2))*(max(a)-min(a)) #R/S statistic


my.loglike.normal=function(theta) #Estimate an asymmetric GARCH(1,1) model
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  q=rep(0,n+1)
  s=rep(0,n+1)
  q[1]=my.sigmasq[1]
  for (i in 2:(n+1))
  {
     q[i]=theta[2]^2*abs(data[i-1]-my.mean[i-1])^2 + theta[3]^2*q[i-1]
     s[i]=theta[4]^2+theta[5]^2*abs(data[i-1]-my.mean[i-1])^2+theta[6]^2*s[i-1]
     my.sigmasq[i]=q[i]+s[i]
  }

  #1/2*sum(log(my.sigmasq[-1])) - sum(log(dnorm((data[-1]-my.mean[-1])/sqrt(my.sigmasq[-1]))))+10^{10}*((theta[2]^2+theta[3]^2)<(theta[5]^2+theta[6]^2))    
   1/2*sum(log(my.sigmasq[2:(n+1)]*(theta[7]-2)/theta[7])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]*(theta[7]-2)/theta[7]),df=theta[7])))+10^(10)*(theta[7]<2)+10^(10)*(theta[7]>10)+10^{10}*((theta[2]^2+theta[3]^2)<(theta[5]^2+theta[6]^2))   
}

returns=ibm
par.start=c(0.5,sqrt(0.03),sqrt(0.9),0.01,sqrt(0.05),sqrt(0.7),4)
my.optpar= nlm(my.loglike.normal,par.start,iterlim=1000,print.level=2)

sigmasq.model=function(theta)
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=sqrt(my.sigmasq)
 
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  q=rep(0,n+1)
  s=rep(0,n+1)
  q[1]=my.sigmasq[1]
  for (i in 2:(n+1))
  {
     q[i]=theta[2]^2*abs(data[i-1]-my.mean[i-1])^2 + theta[3]^2*q[i-1]
     s[i]=theta[4]^2+theta[5]^2*abs(data[i-1]-my.mean[i-1])^2+theta[6]^2*s[i-1]
     my.sigmasq[i]=q[i]+s[i]
  }


	list(my.sigmasq = my.sigmasq[2:(n + 1)],my.mean=my.mean[2:(n+1)])
}

a=sigmasq.model(my.optpar$estimate)

par(mfrow=c(2,2))
acf((ibm-a$my.mean)/sqrt(a$my.sigmasq),main="Residuals",lag.max=50)
acf((ibm-a$my.mean)^2/a$my.sigmasq,main="Squared residuals",lag.max=50)
qqplot((ibm-a$my.mean)/sqrt(a$my.sigmasq), rt(length(ibm), df = my.optpar$estimate[7]),xlab="Theoretical Quantiles",ylab="Sample Quantiles")

b=(ibm-a$my.mean)/sqrt(a$my.sigmasq)
a=c()
for (i in 1:length(ibm))
{
  a=c(a,sum(b[1:i]^2-mean(b^2)))
}

1/sqrt(length(b))*1/sqrt(var(b^2))*(max(a)-min(a)) #R/S statistic

################################################

Chapter 4, page 6ss:

#hk.jp=read.table("d-hkjp0608.txt",header=T)

library(zoo)
a=hk.jp[,3]*10^4+hk.jp[,2]*10^2+hk.jp[,1]
a=strptime(a,"%Y%d%m")
a=a[-1]
hk=log(hk.jp[2:length(hk.jp[,4]),4]/hk.jp[1:(length(hk.jp[,4])-1),4])
jp=log(hk.jp[2:length(hk.jp[,4]),5]/hk.jp[1:(length(hk.jp[,4])-1),5])

par(mfrow=c(2,1))
plot(zoo(hk,a),main="Hang Seng index",ylab="Log returns",xlab="")
plot(zoo(jp,a),main="Nikkei 225 index",ylab="Log returns",xlab="")

library(fGarch) #Load the package for ARCH/GARCH estimation
m1=garchFit(jp~garch(1,1),data=jp,trace=F)
summary(m1) #Obtain results

##EWMA estimation:

r=cbind(hk-mean(hk),jp-mean(jp))
k=dim(r)[2]
n=dim(r)[1]

my.loglike.normal=function(theta)
{   
   residu=t(cbind(t(t(colMeans(r))),t(r)))
   Gamma=var(r)

   l1=0
   l2=0
   for (i in 1:(n))
   {
       Gamma=(1-theta[1]^2)*residu[i,]%*%t(residu[i,])+theta[1]^2*Gamma
	 l1=l1+residu[i+1,]%*%solve(Gamma)%*%residu[i+1,]
	 l2=l2+det(Gamma)
   }

    (n)*k/2*log(2*pi) + 1/2*l2 + 1/2*l1 + 10^(10)*(theta[1]>=1)
}

my.optpar=nlm(my.loglike.normal,0.5,iterlim=1000,print.level=2)

sigmasq.model=function(theta)
{
   residu=t(cbind(t(t(colMeans(r))),t(r)))
   
   V=array(dim=c(n,k,k))
   Gamma=var(r)
  
   for (i in 1:(n))
   {
       V[i,,]=Gamma=(1-theta[1]^2)*residu[i,]%*%t(residu[i,])+theta[1]^2*Gamma
   }

   list(V=V)

}

b=sigmasq.model(my.optpar$estimate)
m1=garchFit(hk~garch(1,1),data=hk,trace=F)
m2=garchFit(jp~garch(1,1),data=jp,trace=F)

par(mfrow=c(2,2)) #Plot of the volatility estimates
plot(zoo(m1@sigma.t,a),main="Hang Seng",xlab="",ylab="GARCH(1,1) volatility")
plot(zoo(m2@sigma.t,a),main="Nikkei 225",xlab="",ylab="GARCH(1,1) volatility")
plot(zoo(sqrt(b$V[,1,1]),a),main="Hang Seng",xlab="",ylab="EWMA volatility")
plot(zoo(sqrt(b$V[,2,2]),a),main="Nikkei 225",xlab="",ylab="EWMA volatility")

#Correlation plot:
plot(zoo(b$V[,1,2]/sqrt(b$V[,1,1]*b$V[,2,2]),a),xlab="",ylab="EWMA Hang Seng / Nikkei correlation")

#Chapter 4, page 12: Diagonal Model

r=cbind(hk-mean(hk),jp-mean(jp))
k=dim(r)[2]
n=dim(r)[1]

my.loglike.normal=function(theta)
{   
   residu=t(cbind(t(t(colMeans(r))),t(r)))
   Gamma=var(r)
 
    A0=matrix(c(3.996e-06,theta[1],theta[1],4.515e-06),k,k)
    A1=matrix(c(1.446e-01,theta[2],theta[2],1.254e-01),k,k)
    B=matrix(c(8.536e-01,theta[3],theta[3],8.629e-01),k,k)

   l1=0
   l2=0
   for (i in 1:(n))
   {
       Gamma=A0+diag(residu[i,])%*%A1%*%diag(residu[i,])+B*Gamma
	 l1=l1+residu[i+1,]%*%solve(Gamma)%*%residu[i+1,]
	 l2=l2+det(Gamma)
   }

    (n)*k/2*log(2*pi) + 1/2*l2 + 1/2*l1 
}

par.start=c(sqrt(3.996e-06*4.515e-06)-10^(-6),sqrt(1.446e-01*1.254e-01)-0.01,sqrt(8.536e-01*8.629e-01)-0.01)
my.optpar=nlminb(par.start,my.loglike.normal,lower=c(0,0,0),upper=c(sqrt(3.996e-06*4.515e-06),sqrt(1.446e-01*1.254e-01),sqrt(8.536e-01*8.629e-01)))

sigmasq.model=function(theta)
{   
   residu=t(cbind(t(t(colMeans(r))),t(r)))
   Gamma=var(r)
 
    A0=matrix(c(3.996e-06,theta[1],theta[1],4.515e-06),k,k)
    A1=matrix(c(1.446e-01,theta[2],theta[2],1.254e-01),k,k)
    B=matrix(c(8.536e-01,theta[3],theta[3],8.629e-01),k,k)

    V=array(dim=c(n,k,k))

   for (i in 1:(n))
   {
       V[i,,]=Gamma=A0+diag(residu[i,])%*%A1%*%diag(residu[i,])+B*Gamma
   }

   list(V=V,A0=A0,A1=A1,B=B)
}

b=sigmasq.model(my.optpar$par)

library(zoo)
a=hk.jp[,3]*10^4+hk.jp[,2]*10^2+hk.jp[,1]
a=strptime(a,"%Y%d%m")
a=a[-1]

plot(zoo(b$V[,1,2]/sqrt(b$V[,1,1]*b$V[,2,2]),a),main="Diagonal model",xlab="",ylab="Correlations")

#################################################

#Chapter 5, page 4:
#source("us.stock.data.R")

par(mfrow=c(2,2))
acf(us.stock.data$returns[,4]^2,lag.max=200,main="Harley Davidson")
acf(us.stock.data$returns[,5]^2,lag.max=200,main="Intel")
acf(us.stock.data$returns[,6]^2,lag.max=200,main="Microsoft")
acf(us.stock.data$returns[,7]^2,lag.max=200,main="Nike")


#Chapter 5, page 9:
par(mfrow=c(2,1))

a=strptime(ibm.sp500[2527:7582,1],"%Y%m%d")
b=sp500[2527:7582]

m2=c() #Monthly real volatility from daily returns
l=c()
u=b[1]^2
v=1
for (i in 2:length(a))
{
   if (months(a)[i-1]==months(a)[i])
    {
        v=v+1
        u=u+b[i]^2
    }
   
   else
    {
        l=c(l,v)
        m2=c(m2,u)
        u=b[i]^2
        v=1
    }

}
l=c(l,v)
m2=c(m2,u)

m=c()  #Check whether monthly returns are compatible with daily returns
u=1
for (i in 1:length(l))
{
  m=c(m,sum(b[u:(u+l[i]-1)]))
  u=u+l[i]
}



b2=ibm.sp500.m[433:888,5] # Monthly returns
#b2=b2-mean(b2)

library(fGarch) 
library(zoo)
m1=garchFit(b2~garch(1,1),data=b2,trace=F)
a=strptime(ibm.sp500.m[433:888,1],"%Y%m%d")

plot(zoo(sqrt(m2),a[217:456]),xlab="Year",ylab="Volatility")
plot(zoo(m1@sigma.t[217:456],a[217:456]),xlab="Year",ylab="Volatility")

#Chapter 5, page 15:

boeing1=read.table("taq-td-ba12012008.txt",header=T) #Load the tick-data for Boeing for the first five days in December 2008
boeing2=read.table("taq-td-ba12022008.txt",header=T)
boeing3=read.table("taq-td-ba12032008.txt",header=T)
boeing4=read.table("taq-td-ba12042008.txt",header=T)
boeing5=read.table("taq-td-ba12052008.txt",header=T)

"hfrtn" = function(da,int,logrtn=TRUE){
# Compute intraday returns
#
# int: time intervals in minutes
# da: data in the format: date, hour, minute, second, price, volume
#
if(!is.matrix(da))da=as.matrix(da)
intsec=int*60
istart=9*60*60+30*60
iend=16*60*60
# compute the number of prices
tradetime=6.5*60*60
ntrade=floor(tradetime/intsec)
T=dim(da)[1]
nday=da[T,1]-da[1,1]+1
npri=nday*ntrade
#print(c(ntrade,nday,npri))

price=rep(0,npri)
# price is the last transaction price of the time interval
caltime=da[,2]*60*60+da[,3]*60+da[,4]
#plot(caltime,type='l')

icnt=0
date=da[1,1]
for (i in 1:T) {
if(caltime[i] &gt; istart){
iday=da[i,1]-date
if(caltime[i] < (iend+1)){

if(caltime[i]==iend){
price[iday*ntrade+ntrade]=da[i,5]
}

if((caltime[i] &gt; istart) &amp;&amp; (caltime[i] < iend)){
ii=caltime[i]-istart
ij=floor(ii/intsec)
price[iday*ntrade+ij+1]=da[i,5]
}
}
}
}
for (i in 2:npri){
if(price[i] <= 0)price[i]=price[i-1]
}

plot(price,type='l')

pri=log(price)
#skip overnight returns
nrtn=ntrade-1
rtn=NULL
for (i in 1:nday){
ist=(i-1)*ntrade
for (j in 2:ntrade){
rtn=c(rtn,pri[ist+j]-pri[ist+j-1])
}
}

hfrtn = list(rtn=rtn,price=price)
}


rv.boeing=matrix(,5,60)

for (t in 1:60)
{
   a=hfrtn(da=boeing5,int=t)

for (i in 1:length(a$rtn))
{
  if (a$rtn[i]=="NaN")
   {
      a$rtn[i]=0
   }
}

rv.boeing[5,t]=sum(a$rtn^2)
}

rv=rv.boeing
a=colMeans(rv)

ts.plot(a,xlab="Minutes",ylab="Sample RV",main="Volatility signature plot") #Volatility signature plot
lines(rep(mean(a[30:60]),60),lty=3)

#Chapter 5, page 18ss:

par(mfrow=c(2,2))
acf(log(us.stock.data$real.cov[,4,4]),lag.max=200,main="Harley Davidson")
acf(log(us.stock.data$real.cov[,5,5]),lag.max=200,main="Intel")
acf(log(us.stock.data$real.cov[,6,6]),lag.max=200,main="Microsoft")
acf(log(us.stock.data$real.cov[,7,7]),lag.max=200,main="Nike")


a=log(us.stock.data$real.cov[,7,7]) #Estimation of HAR model for log(RV):

u=21 #Number of past lags needed in HAR model
n=length(a)-u
dat.pred=matrix(,n,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(length(a)-u+i-1)]
} 
dat.resp=a[(u+1):length(a)]

v=5 #Horizon for the prediction

returns=c()
for (i in v:n)
{
  returns=c(returns,mean(dat.resp[(i-v+1):i]))
}

predictors=matrix(,length(returns),3)
predictors[,1]=dat.pred[1:(n-v+1),1]
predictors[,2]=rowMeans(dat.pred[1:(n-v+1),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n-v+1),1:21])

ff=lm(returns~predictors)
summary(ff)

#############################################

#Chapter 6, page 8:

a=us.stock.data$real.cov[,6,6] #Daily Microsoft RV proxies from 10-minute returns

u=21 #Number of past lags needed in HAR model
n=length(a)-u
dat.pred=matrix(,n,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(length(a)-u+i-1)]
} 
dat.resp=a[(u+1):length(a)]

returns=us.stock.data$returns[(u+1):length(a),6]

my.loglike.normal=function(theta) #Estimate an (asymmetric) GARCH(1,1) model
{
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,length(returns)+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,length(returns)))
  log.sigmasq=c(log(my.sigmasq[1]),rep(0,length(returns)))
  
  my.mean=rep(0,length(returns)+1)
  for(j in 2:(length(returns)+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(length(returns)+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    #aiuto=(data[i-1]-my.mean[i-1])/sqrt(exp(log.sigmasq[i-1]))
    #log.sigmasq[i]=theta[2]+theta[3]*(theta[4]*aiuto+abs(aiuto)-sqrt(2/pi))+theta[5]*log.sigmasq[i-1] #EGARCH(1,1) 
    
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)
  1/2*sum(log(my.sigmasq[-1])) - sum(log(dnorm((data[-1]-my.mean[-1])/sqrt(my.sigmasq[-1]))))
 
}

par.start=rep(0.5,5)
my.optpar= nlm(my.loglike.normal,par.start,iterlim=1000,print.level=2)


sigmasq.model=function(theta)
{
  x.start= mean(returns)
  sigmasq.start= var(returns)
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,length(returns)+1)
  my.sigmasq[1]=sigmasq.start

  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,length(returns)))
  log.sigmasq=c(log(my.sigmasq[1]),rep(0,length(returns)))
  
  my.mean=rep(0,length(returns)+1)
  for(j in 2:(length(returns)+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(length(returns)+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)

    #aiuto=(data[i-1]-my.mean[i-1])/sqrt(exp(log.sigmasq[i-1]))
    #log.sigmasq[i]=theta[2]+theta[3]*(theta[4]*aiuto+abs(aiuto)-sqrt(2/pi))+theta[5]*log.sigmasq[i-1] #EGARCH(1,1) 
    
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2*(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }

  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)


	list(my.sigmasq = my.sigmasq[-1],my.mean=my.mean[-1])
}

s=sigmasq.model(my.optpar$estimate)

rv.proxy=dat.resp
#rv.proxy=returns^2

par(mfrow=c(1,3))
plot(sqrt(s$my.sigmasq),sqrt(rv.proxy),xlab="GARCH(1,1) predictions",ylab="Volatility proxy")
abline(0,1,lty=3)

plot(sqrt(v$my.sigmasq),sqrt(rv.proxy),xlab="GJR-GARCH(1,1) predictions",ylab="Volatility proxy")
abline(0,1,lty=3)

library(car)

b=lm(rv.proxy~v$my.sigmasq) #Univariate in-sample MZ test:
summary(b)
vcov(b) #OLE standard errors
hccm(b,type="hc0") #White consistent standard errors

b=lm((rv.proxy/v$my.sigmasq)[2:n]~(rv.proxy/v$my.sigmasq)[1:(n-1)]) #Alternative test 
summary(b)


predictors=matrix(,length(dat.resp),3) #HAR estimation
predictors[,1]=dat.pred[1:(n),1]
predictors[,2]=rowMeans(dat.pred[1:(n),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n),1:21])

m1=lm(dat.resp~predictors) 

plot(sqrt(m1$fitted.values),sqrt(rv.proxy),xlab="HAR predictions",ylab="Volatility proxy")
abline(0,1,lty=3)


b=lm(rv.proxy~m1$fitted.values) #Univariate in-sample MZ test:
summary(b)
vcov(b) #OLE standard errors
hccm(b,type="hc0") #White consistent standard errors

b=lm((rv.proxy/m1$fitted.values)[2:n]~(rv.proxy/m1$fitted.values)[1:(n-1)]) #Alternative test 
summary(b)

par(mfrow=c(3,1))
plot(zoo(sqrt(dat.resp),us.stock.data$dates[-c(1:21)]),main="Realized volatility",xlab="Year",ylab="")
plot(zoo(sqrt(returns^2),us.stock.data$dates[-c(1:21)]),main="Absolute daily returns",xlab="Year",ylab="")
plot(zoo(sqrt(s$my.sigmasq),us.stock.data$dates[-c(1:21)]),main="GARCH(1,1) volatility",xlab="Year",ylab="")


#Chapter 6, page 12:

library(car)

qt(0.99,df=2462-1) #p-value

for (s in 4:7)
{
for (t in s:7)
{
print(s)
print(t)

a=us.stock.data$real.cov[,s,t] 

u=21 #Number of past lags needed in HAR model
n=length(a)-u
dat.pred=matrix(,n,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(length(a)-u+i-1)]
} 
dat.resp=a[(u+1):length(a)]

returns=us.stock.data$returns[(u+1):length(a),6]


v=us.stock.data$dcc.estimates[-c(1:u),s,t] #Previously estimated CCC/DCC-GARCH(1,1) covariances

rv.proxy=dat.resp
#rv.proxy=us.stock.data$returns[(u+1):length(a),s]*us.stock.data$returns[(u+1):length(a),t]

b=lm(rv.proxy~v) #Univariate in-sample MZ test:
print(b$coef[1]/sqrt(hccm(b,type="hc0")[1,1]))
print((1-b$coef[2])/sqrt(hccm(b,type="hc0")[2,2]))


predictors=matrix(,length(dat.resp),3) #HAR estimation
predictors[,1]=dat.pred[1:(n),1]
predictors[,2]=rowMeans(dat.pred[1:(n),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n),1:21])

m1=lm(dat.resp~predictors) 

b=lm(rv.proxy~m1$fitted.values) #Univariate in-sample MZ test:
print(b$coef[1]/sqrt(hccm(b,type="hc0")[1,1]))
print((1-b$coef[2])/sqrt(hccm(b,type="hc0")[2,2]))

}
}

#Joint MZ test:

library(systemfit)

har=matrix(,n,10)
garch=matrix(,n,10)
rv.proxy=matrix(,n,10)

v=1
for (s in 4:7)
{
for (t in s:7)
{
a=us.stock.data$real.cov[,s,t] 

u=21 #Number of past lags needed in HAR model
n=length(a)-u
dat.pred=matrix(,n,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(length(a)-u+i-1)]
} 
dat.resp=a[(u+1):length(a)]

predictors=matrix(,length(dat.resp),3) #HAR estimation
predictors[,1]=dat.pred[1:(n),1]
predictors[,2]=rowMeans(dat.pred[1:(n),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n),1:21])

har[,v]=lm(dat.resp~predictors)$fitted.values 
garch[,v]=us.stock.data$dcc.estimates[-c(1:u),s,t]
#rv.proxy[,v]=dat.resp
rv.proxy[,v]=us.stock.data$returns[(u+1):length(a),s]*us.stock.data$returns[(u+1):length(a),t]

v=v+1
}
}

system=list(eq1=rv.proxy[,1]~garch[,1],eq2=rv.proxy[,2]~garch[,2],eq3=rv.proxy[,3]~garch[,3],eq4=rv.proxy[,4]~garch[,4],eq5=rv.proxy[,5]~garch[,5],eq6=rv.proxy[,6]~garch[,6],eq7=rv.proxy[,7]~garch[,7],eq8=rv.proxy[,8]~garch[,8],eq9=rv.proxy[,9]~garch[,9],eq10=rv.proxy[,10]~garch[,10])

m1=systemfit(system)
linearHypothesis(m1,hypothesis.matrix=diag(20),rhs=rep(c(0,1),10),test = "Chisq")


#Chapter 6, page 15:

library(mvtnorm)

mod=us.stock.data$ccc.estimates

a=0  #QLIKE
for (i in 1:length(mod[,1,1]))
{
 a=a-log(dmvnorm(us.stock.data$returns[i,4:7],mean=colMeans(us.stock.data$returns[,4:7]),sigma=mod[i,4:7,4:7]))
}
print(a)

a=matrix(,4,4)  #MAE/MSE
b=matrix(,4,4)

for (i in 4:7)
{
	for (j in 4:7)
	{   
		a[i-3,j-3]=mean(abs(us.stock.data$real.cov[,i,j]-mod[,i,j]))
		a[j-3,i-3]=a[i-3,j-3]
		b[i-3,j-3]=mean((us.stock.data$real.cov[,i,j]-mod[,i,j])^2)
        b[j-3,i-3]=b[i-3,j-3]
	}
}

mean(a)
mean(b)

#DMW test:

n=length(mod[,1,1])
perf=rep(0,length(mod[,1,1])) 
for (i in 1:length(mod[,1,1]))
{
perf[i]= length(mod[,1,1])*(-log(dmvnorm(us.stock.data$returns[i,4:7],mean=colMeans(us.stock.data$returns[,4:7]),sigma=us.stock.data$ccc.estimates[i,4:7,4:7]))+log(dmvnorm(us.stock.data$returns[i,4:7],mean=colMeans(us.stock.data$returns[,4:7]),sigma=us.stock.data$dcc.estimates[i,4:7,4:7]))) #QLIKE
#perf[i]= mean(abs(us.stock.data$real.cov[i,4:7,4:7]-us.stock.data$ccc.estimates[i,4:7,4:7])) - mean(abs(us.stock.data$real.cov[i,4:7,4:7]-us.stock.data$dcc.estimates[i,4:7,4:7])) #MAE
#perf[i]= mean(abs(us.stock.data$real.cov[i,4:7,4:7]-us.stock.data$ccc.estimates[i,4:7,4:7])^2) - mean(abs(us.stock.data$real.cov[i,4:7,4:7]-us.stock.data$dcc.estimates[i,4:7,4:7])^2) #MSE
}

sqrt(n)*mean(perf)/sqrt(spectrum(perf)$spec[1]) 
1-pnorm(sqrt(n)*mean(perf)/sqrt(spectrum(perf)$spec[1]))


#Chapter 6, page 17:

#Create u a vector with all different univariate GARCH predictions: each column, one model prediction

n=dim(u)[1]
perf=matrix(,n,5)

for (j in 1:5)
{
  #perf[,j]= n*(-log(dnorm(us.stock.data$returns[,6],mean=mean(us.stock.data$returns[,6]),sd=sqrt(u[,j])))) #QLIKE
  #perf[,j]= abs(us.stock.data$real.cov[,6,6]-u[,j]) #MAE
  perf[,j]= abs(us.stock.data$real.cov[,6,6]-u[,j])^2 #MSE
}

colMeans(perf)

spa=function(per=perf,bench=1,m=9,obs=926,q=0.25,iter=1,periodogram=T) #SPA test
{
#Test of superior predictive ability of Hansen, 2005 JBES
e=bench #benchmark
d=matrix(,obs,m-1)
s=0
for (i in seq(1,m,1)[-e])
{
	s=s+1
	d[,s]=per[,e]-per[,i]
}

#colMeans(d)

w=rep(0,m-1)
for (k in 1:(m-1))
{
  #e=c()
  #for (i in 0:(obs-1))
  #{
     #e=c(e,1/obs*sum((d[1:(obs-i),k]-mean(d[,k]))*(d[(1+i):(obs),k]-mean(d[,k]))))
      e=acf(d[,k],lag.max=obs-1,type="covariance",plot=F)$acf
  #}
  if (periodogram==F)
  {
	w[k]=sqrt(e[1]+2*sum(((obs-seq(1,obs-1,1))/obs*(1-q)^{seq(1,obs-1,1)}+seq(1,obs-1,1)/obs*(1-q)^{obs-seq(1,obs-1,1)})*e[2:obs]))
  }
  else if (periodogram==T)
 {
  w[k]=sqrt(spectrum(d[,k],plot=F)$spec[1])
 }
}

#print(sqrt(obs)*colMeans(d)/w)
#print(pnorm(sqrt(obs)*colMeans(d)/w))

stat=max(0,max(sqrt(obs)*colMeans(d)/w))

#Bootstrap:

stat.boos=rep(0,iter)
for (r in 1:iter)
{
#print(r)
tau=rep(0,obs)
tau[1]=as.integer(obs*runif(1))+1
for (i in 2:obs)
{
	s=runif(1)
	tau[i]=(as.integer(obs*runif(1))+1)*(s<q)+((tau[i-1]<obs)*tau[i-1]+1)*(s&gt;=q)
}

d.boos=d[tau,]

e=d
for (k in 1:(m-1))
{
  e[,k]=d.boos[,k]-mean(d[,k])*(mean(d[,k])&gt;= - sqrt(w[k]^2/obs*2*log(log(obs))))
}

stat.boos[r]=max(0,max(sqrt(obs)*colMeans(e)/w))
}
p.value=mean((stat.boos&gt;stat))

list(p.value=p.value,stat.boos=stat.boos,stat=stat)
}


for (k in 1:5)
{
print(k)	
d=rep(0,1000)
s=1
for (i in 1:10)
{ 
	#print(i)
   e=spa(per=perf,bench=k,m=5,obs=n,q=0.25,iter=100,periodogram=T)
   d[s:(s+99)]=e$stat.boos
   s=s+100
}
print(mean((d&gt;e$stat)))
}


#Chapter 6, page 23

library(quadprog)

n1=1507
n2=976
k=9 #Number of asssets

a=matrix(,n2,k)
b=matrix(,n2,k)
d=matrix(,n2,k)

perf=matrix(,n2,3)

for (t in (n1+1):(n1+n2))
{
   a[t-n1,] = solve.QP(Dmat=us.stock.data$ccc.estimates[t,,], dvec=array(0, dim = c(k,1)), Amat=t(array(1, dim = c(1,k))), bvec=1, meq = 1)$solution #Global minimum variance portfolio
   b[t-n1,] = solve.QP(Dmat=us.stock.data$dcc.estimates[t,,], dvec=array(0, dim = c(k,1)), Amat=t(array(1, dim = c(1,k))), bvec=1, meq = 1)$solution
   d[t-n1,] = solve.QP(Dmat=us.stock.data$real.cov[t-1,,], dvec=array(0, dim = c(k,1)), Amat=t(array(1, dim = c(1,k))), bvec=1, meq = 1)$solution
   
   #perf[t-n1,1]=t(a[t-n1,])%*%t(us.stock.data$returns[t,])%*%t(t(us.stock.data$returns[t,]))%*%a[t-n1,] #With cross returns as proxies
   #perf[t-n1,2]=t(b[t-n1,])%*%t(us.stock.data$returns[t,])%*%t(t(us.stock.data$returns[t,]))%*%b[t-n1,]
   #perf[t-n1,3]=t(d[t-n1,])%*%t(us.stock.data$returns[t,])%*%t(t(us.stock.data$returns[t,]))%*%d[t-n1,]

   perf[t-n1,1]=t(a[t-n1,])%*%us.stock.data$real.cov[t,,]%*%a[t-n1,] #With real vol as proxies
   perf[t-n1,2]=t(b[t-n1,])%*%us.stock.data$real.cov[t,,]%*%b[t-n1,]
   perf[t-n1,3]=t(d[t-n1,])%*%us.stock.data$real.cov[t,,]%*%d[t-n1,]

}

colMeans(perf)

perf=perf[,2]-perf[,3]

sqrt(n2)*mean(perf)/sqrt(spectrum(perf)$spec[1])  #DMW test
1-pnorm(sqrt(n2)*mean(perf)/sqrt(spectrum(perf)$spec[1]))


#Chapter 6, page 30:
library(VGAM)
a=c()
a[41:1]=seq(0.95,0.99,0.001)

plot(seq(0.01,0.05,0.001),qnorm(a),type="l",ylim=c(1.4,2.8),xlab="alpha",ylab="VaR")
lines(seq(0.01,0.05,0.001),sqrt(1/3)*qt(a,df=3),lty=2)
lines(seq(0.01,0.05,0.001),qlaplace(a, location=0, scale=sqrt(1/2)),lty=3)


#Chapter 6, page 39:

returns=us.stock.data$returns[,4]/100 -mean(us.stock.data$returns[,4]/100)
b=100
prices=b
for (i in 1:n2)
{
  b=b*exp(returns[n1+i])
  prices=c(prices,b)
}

alpha=0.01

#Unconditional VaR:

a=c()
for (i in 1:n2)
{
   a=c(a,(1-exp(quantile(returns[(n1+i-1):(n1+i-1-250+1)],alpha)))*prices[i])    
}

mean(a)

plot(-(prices[2:(n2+1)]-prices[1:n2]),type="l",main="Unconditional VaR",xlab="",ylab="",ylim=c(-9,9))
lines(a,lty=3)

sum(-(prices[2:(n2+1)]-prices[1:n2])&gt; a) #Number of violations

#RiskMetrics:

b=var(returns[1:n1])
for (i in 1:n2)
{ 
   b=c(b,0.94*b[i]+(1-0.94)*returns[n1+i-1]^2)
}

a=(1-exp(sqrt(b[-1])*qnorm(alpha)))*prices[1:n2]

mean(a)

plot(-(prices[2:(n2+1)]-prices[1:n2]),type="l",main="RiskMetrics VaR",xlab="",ylab="",ylim=c(-9,9))
lines(a,lty=3)

sum(-(prices[2:(n2+1)]-prices[1:n2])&gt; a) #Number of violations

#GARCH

library(fGarch) #Load the package for ARCH/GARCH estimation
m1=garchFit(returns[1:n1]~garch(1,1),data=returns[1:n1],trace=F)

sigmasq.model=function(theta)
{
  n=length(returns)
  x.start= mean(returns[1:n1])
  sigmasq.start= var(returns[1:n1])
   
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start

  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  
  for (i in 2:(n+1))
  {
    my.sigmasq[i]=theta[2] + theta[3]*(data[i-1]-my.mean[i-1])^2 + theta[4]*my.sigmasq[i-1] #GARCH(1,1)
  }

	list(my.sigmasq = my.sigmasq[2:(n + 1)],my.mean=my.mean[2:(n+1)])
}

m1=sigmasq.model(m1@fit$matcoef[,1])

a=(1-exp(sqrt(m1$my.sigmasq[(n1+1):(n1+n2)])*qnorm(alpha)))*prices[1:n2]  #Under a Gaussian distribution

mean(a)

plot(-(prices[2:(n2+1)]-prices[1:n2]),type="l",main="GARCH(1,1)-N VaR",xlab="",ylab="",ylim=c(-9,9))
lines(a,lty=3)

sum(-(prices[2:(n2+1)]-prices[1:n2])&gt; a) #Number of violations


my.loglike.t=function(theta) 
{

  x.start= mean(returns[1:n1])
  sigmasq.start= var(returns[1:n1])
   
  data=c(x.start,returns[1:n1])
  my.sigmasq= rep(0,n1+1)
  my.sigmasq[1]=sigmasq.start

  my.mean=rep(0,n1+1)
  for(j in 2:(n1+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }

  for (i in 2:(n1+1))
  {
   my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)
  }

  1/2*sum(log(my.sigmasq[2:(n1+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n1+1)]-my.mean[2:(n1+1)])/sqrt(my.sigmasq[2:(n1+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]&gt;10)
 }

par.start=c(rep(0.5,5),3)
my.optpar= nlm(my.loglike.t,par.start,iterlim=1000,print.level=2)

m1=sigmasq.model(c(my.optpar$estimate[1],my.optpar$estimate[2:4]^2))

a=(1-exp(m1$my.mean[(n1+1):(n1+n2)]+sqrt(m1$my.sigmasq[(n1+1):(n1+n2)])*qt(alpha,my.optpar$estimate[6])))*prices[1:n2]  #Under a Students t-distribution

mean(a)

plot(-(prices[2:(n2+1)]-prices[1:n2]),type="l",main="GARCH(1,1)-t VaR",xlab="",ylab="",ylim=c(-9,9))
lines(a,lty=3)

sum(-(prices[2:(n2+1)]-prices[1:n2])&gt; a) #Number of violations

#HAR:

a=log(us.stock.data$real.cov[,4,4]/100^2) 

u=21 #Number of past lags needed in HAR model
dat.pred=matrix(,n1-u,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(n1-u+i-1)]
} 
dat.resp=a[(u+1):n1]

predictors=matrix(,length(dat.resp),3) 
predictors[,1]=dat.pred[1:(n1-u),1]
predictors[,2]=rowMeans(dat.pred[1:(n1-u),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n1-u),1:21])

b=lm(dat.resp~predictors) #HAR estimation in-sample

dat.pred=matrix(,n-u,u)
for (i in 1:u)
{
   dat.pred[,u-i+1]=a[i:(n-u+i-1)]
} 
dat.resp=a[(u+1):n]

predictors=matrix(,length(dat.resp),3) #HAR estimation
predictors[,1]=dat.pred[1:(n-u),1]
predictors[,2]=rowMeans(dat.pred[1:(n-u),1:5])
predictors[,3]=rowMeans(dat.pred[1:(n-u),1:21])

m1=c()
for (i in 1:n2)
{
  m1=c(m1,b$coef[1]+b$coef[2]*predictors[n1-u+i,1]+b$coef[3]*predictors[n1-u+i,2]+b$coef[4]*predictors[n1-u+i,3])
}


a=(1-exp(sqrt(exp(m1))*qnorm(alpha)))*prices[1:n2]  #Under a Gaussian distribution

mean(a)

plot(-(prices[2:(n2+1)]-prices[1:n2]),type="l",main="HAR VaR",xlab="",ylab="",ylim=c(-9,9))
lines(a,lty=3)

sum(-(prices[2:(n2+1)]-prices[1:n2])&gt; a) #Number of violations













