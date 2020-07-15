setwd("C:/Users/charr/Desktop/RMIT/#MATH1318 Time Series Analysis/A3 - May 29")

# INSTALLATION AND LOADING OF THE TIME SERIES PACKAGES
ts_packages <- c('TSA', 'fUnitRoots', 'forecast', 'lmtest', 'FitAR', 'tseries')
for (pkg in ts_packages) {
  if (pkg %in% rownames(installed.packages()) == FALSE)
  {install.packages(pkg)}
  if (pkg %in% rownames(.packages()) == FALSE)
  {library(pkg, character.only = TRUE)}
}

# FUNCTION-SETTING

# Descriptive Analysis Function
descriptive_analysis <- function(ts, ylabel="ylabel", xlabel="xlabel", title="Name of the Time Series") {
  # Time Series Plot
  plot(ts, ylab=ylabel, xlab=xlabel, type="o", main=paste("Time Series Plot of", title))
  # Correlation Analysis
  plot(y=ts,x=zlag(ts), ylab=ylabel, xlab = "Previous Month Change", main = paste("Scatter Plot of ", title))
  y = ts
  x = zlag(ts) 
  index = 2:length(x)
  cor(y[index],x[index])
}

# Normality Test function
normality_test <- function(x, title = "Series Name") { 
  par(mfrow=c(1,2))
  qqnorm(x, main="")
  qqline(x, col = "blue")
  data_sd <- sd (x)
  data_mean <- mean(x)
  data_median <- median(x)
  hist (x, main = "", probability = TRUE, xlab = "Time Series Data")
  xm <- seq(min(x),max(x),length = 100)
  ym <- dnorm (x= xm, mean = data_mean, sd = data_sd)
  lines(x = xm, y= ym, col = "blue", lwd = 1.5)
  abline (v = data_mean, col = "red", lwd = 3, lty = 2)
  abline (v = data_median, col = "dark green", lwd = 3, lty = 2)
  legend("topright", legend = c("Mean","Median"), col = c("red", "dark green"), pch = 15, bty = "n")
  par(mfrow=c(1,1))
  mtext(paste("Normal Q-Q Plot and Histogram of", title), side = 3, line = -2, outer = TRUE)
  print(shapiro.test(x))
}

# Deterministic trend function
deterministic_model <- function(x, model = c("Linear", "Quadratic"), title="Series Name", ylabel="Label", xlabel="Label") {
  t = time(x)
  t2 = t^2
  
  if (model == "Linear") {
    model_x = lm(x~t)
    summary_model <- summary(model_x)
    plot_model <- { plot(x, type='o', ylab=ylabel, xlab=xlabel, main = paste(model, "Trend Model of", title))
      abline(model_x, col ="red")
      legend("topright", lty =1, bty = "n", col = c("black", "red"), c("Time Series Plot", "Model fitted value"))}
    res_analysis <- {
      # Standardized Residuals vs. Time
      res.model_x = rstudent(model_x)
      plot(y = res.model_x, x = as.vector(time(x)),xlab=xlabel, ylab='Standardized Residuals', type='p', main="Scatter Plot of the Standarized Residuals")
      abline(h=0)
      # Standardised residuals vs. fitted trend values
      scatter.smooth(y=res.model_x, x=fitted(model_x), xlab='Fitted Trend Values', ylab='Standardized Residuals', main = "Plot of the Standardised Residuals vs. Fitted Trend Values")
      abline(h = 0, lty = 2)
      Acf(res.model_x, main="ACF of the Standardized Residuals")
      # Test constant variances
      print(ncvTest(model_x))
      # Test residual independence
      print(runs(res.model_x))
      # Normality Test
      normality_test(res.model_x, title) }
  } else if (model == "Quadratic") {
    model_x = lm(x~ t + t2)
    summary_model <- summary(model_x)
    plot_model <- { plot(ts(fitted(model_x)), ylim = c(min(c(fitted(model_x),as.vector(x))), max(c(fitted(model_x),as.vector(x)))), col="red", ylab=ylabel, xlab=xlabel, main = paste(model, " Trend Model of", title))
      lines(as.vector(x), type="o")
      legend("bottomleft", lty =1, bty = "n", text.width = 12, col = c("black", "red"), c("Time Series plot", "Model fitted values"))  
    res_analysis <- {
      # Standardized Residuals vs. Time
      res.model_x = rstudent(model_x)
      plot(y = res.model_x, x = as.vector(time(x)),xlab =xlabel, ylab='Standardized Residuals', type='p', main="Scatter Plot of the Standarized Residuals")
      abline(h=0)
      # Standardised residuals vs. fitted trend values
      scatter.smooth(y=res.model_x,x=fitted(model_x), xlab='Fitted Trend Values', ylab='Standardized Residuals', main = "Plot of the Standardized Residuals vs. Fitted Trend Values")
      abline(h = 0, lty = 2)
      Acf(res.model_x, main="ACF of the Standardized Residuals")
      # Test constant variances
      print(ncvTest(model_x))
      # Test residual independence
      print(runs(res.model_x))
      # Normality Test
      normality_test(res.model_x, title) }
    }
  } else { stop("Choose between Linear or Quadratic only")
  }
  return(summary_model)
  return(plot_model)
  return(res_analysis)
  par(mfrow=c(1,1))
}

# ACF and PACF plot function
correlograms <- function(x, title = "Series Name") {
  par(mfrow=c(1,3))
  Acf(x, main = "")
  acf(x, ci.type='ma', main= "")
  Pacf(x, main = "")
  mtext(paste("ACF/PACF Plot of", title), side = 3, line = -2, outer = TRUE)
  par(mfrow=c(1,1))
}

# Differencing function
differencing <- function(x){
  # First Differencing
  data <- diff(x, diff = 1)
  order <-  ar(diff(data))$order
  p <- adfTest(data, lags = order,  title = NULL,description = NULL)@test$p.value
  
  # Second Differencing
  if (p >= 0.05) {
    data_2df <- diff(x, diff = 2)
    order <- ar(diff(data_2df))$order
    p2 <- adfTest(data_2df, lags = order,  title = NULL,description = NULL)@test$p.value
    p2
    if (p2 >= 0.05) {print("Model may result to overdiffercing after Second Differencing. Check again.")
    }else {
      print(shapiro.test(data_2df))
      print("Achieved stationarity after Second Differencing")
      print(paste("p-value from the adfTest", p2))
    }
  }else if (p < 0.05) {
    print(shapiro.test(data))
    print("Achieved stationarity after First Differencing")
    print(paste("p-value from the adfTest", p))
  }else {stop("Check error")
  }
}

# Sort score function
sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

# Residual analysis function
residual.analysis <- function(model, std = TRUE){
  if(std == TRUE){
    res.model = rstandard(model)
  }else{
    res.model = residuals(model)
  }
  par(mfrow=c(3,2))
  plot(res.model, type = 'o', ylab = 'Standardised residuals', main = "Time Series Plot of the Standardised Residuals")
  abline(h = 0)
  hist(res.model, main = "Histogram of the Standardised Residuals")
  qqnorm(res.model, main = "Normal Q-Q Plot of the Standardised Residuals")
  qqline(res.model, col = 2)
  Acf(res.model, main = "ACF of standardised residuals")
  print(shapiro.test(res.model))
  print(runs(res.model))
  k=0
  LBQPlot(res.model, lag.max = length(model$residuals)-1, StartLag = k+1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
  boxtest <- Box.test(res.model, lag = 12, type = "Ljung-Box", fitdf = 0) # Lags limited to 12 as per the highest df in the sort score
  print(boxtest)
}

# AICc function
AICc = function(model){
  n = model$nobs
  k = length(model$coef)
  aicc = model$aic + 2*(k+1)*(k+2)/(n-k-2)
  return(aicc)
}

# ---------------------------
# READ DATA AND CONVERT TO TIME SERIES DATA
library(readr)
library(tidyr)
library(car)
USunemployment <- read_csv("USUnemployment.csv")
head(USunemployment)
unemp_filter <- subset(USunemployment, Year>=2000)
USunem <- gather(unemp_filter, "month", "unemployment rate", na.rm = TRUE, 2:13)
head(USunem)
USunem_sorted <- USunem[order(USunem$Year),]
head(USunem_sorted)
unem <- ts(USunem_sorted$`unemployment rate`, start=c(2000,1), end=c(2019,12), frequency = 12)

# CHAPTER 1: MONTHLY U.S. UNEMPLOYMENT RATE FROM 2000-2019

# 1.1. DESCRIPTIVE ANALYSIS
descriptive_analysis(unem, "Unemployment Rate (%)", "Year", "the Monthly U.S. Unemployment Rate for 2000-2019")
Acf(unem, main = "ACF of the Monthly U.S. Unemployment Rate from 2000 to 2019")

# 1.2. LINEAR TREND MODEL
deterministic_model(unem, model="Linear", title="the Monthly U.S. Unemployment Rate (2000-2019)", ylabel="Unemployment Rate (%)", xlabel="Year")

# 1.3. QUADRATIC TREND MODEL
par(mfrow=c(1,1))
deterministic_model(unem, model="Quadratic", title="the Monthly U.S. Unemployment Rate (2000-2019)", ylabel="Unemployment Rate (%)", xlabel="Year")

# 1.4. STOCHASTIC TREND MODEL

# 1.4.1 Test of Normality and Stationarity
normality_test(unem, "the U.S. Unemployment Rate")
correlograms(unem, "the U.S. Unemployment Rate")
adf.test(unem)
BoxCox(unem, interval = c(-1, 1)) # Data Transformation
lambda = -0.709
unem.BC = (unem^lambda-1)/lambda
normality_test(unem.BC,  "the Time Series (2000-2019) after Box-Cox Transformation")
order = ar(diff(unem.BC))$order
adfTest(unem.BC, lags = order, type = 'nc',  title = NULL, description = NULL)
plot(unem.BC, ylab='Unemployment Rate (%)',xlab='Year',type='o', main = "Time Series Plot of the Monthly U.S. Unemployment Rate after Box-Cox Transformation")

# 1.4.2. Model Specification
differencing(unem.BC)
diff.unem.BC = diff(unem.BC, diff = 1) # Set data to first differencing
plot(diff.unem.BC ,ylab='Unemployment Rate (%)',xlab='Year',type='o', main = "Time Series Plot of the First Differenced Series (2000-2019)")
abline(h=0, col = 'red', lty = 2)

correlograms(diff.unem.BC, "the First Differenced Series (2000-2019)")
eacf(diff.unem.BC, ar.max = 10, ma.max= 10)
bic = armasubsets(y=diff.unem.BC,nar=10,nma=10,y.name='test',ar.method='ols')
plot(bic)

# 1.4.3. Model Fitting and Parameter Estimation
model.115 <-  arima(unem.BC, order=c(1,1,5), method='ML')
coeftest(model.115)
model.516 <-  arima(unem.BC, order=c(5,1,6), method='ML')
coeftest(model.516)
model.315 <- arima(unem.BC, order=c(3,1,5), method='ML')
coeftest(model.315)
model.316 <- arima(unem.BC, order=c(3,1,6), method='ML')
coeftest(model.316)
model.415 <- arima(unem.BC, order=c(4,1,5), method='ML')
coeftest(model.415)
model.710 <- arima(unem.BC, order=c(7,1,0), method='ML')
coeftest(model.710)
model.711 <- arima(unem.BC, order=c(7,1,1), method='ML')
coeftest(model.711)

# 1.4.4. Residual Analysis
residual.analysis(model.115)
residual.analysis(model.516)
residual.analysis(model.315)
residual.analysis(model.316)
residual.analysis(model.415)
residual.analysis(model.710)
residual.analysis(model.711)

# 1.4.5. Model Selection (of the remaining candidate models)
sort.score(AIC(model.115, model.415), score = "aic")
sort.score(BIC(model.115, model.415), score = "bic")

AICc(model.115)
AICc(model.415)

# 1.4.6. Overfitting (of the best model)
model.215 = arima(unem.BC, order=c(2,1,5), method='ML')
coeftest(model.215)
model.116 = arima(unem.BC, order=c(1,1,6), method='ML')
coeftest(model.116)

# 1.5. FORECAST
fit_115 = Arima(unem,c(1,1,5), lambda = -0.709)
predict_115 = forecast(fit_115, h=10)
predict_115
plot(predict_115,  type = 'o', colour = "blue", xlab = "Year", ylab = "Unemployment Rate (%)", main = "Time Series Plot of the U.S. Unemployment Rate with a Ten-month Forecast") 

# ---------------------------
# CHAPTER 2: MONTHLY U.S. UNEMPLOYMENT RATE FROM 2010-2019

# DATA PREPROCESSING
unemp2 <- subset(USunemployment, Year>=2010)
unem2_gather <- gather(unemp2, "month", "unemployment rate", na.rm = TRUE, 2:13)
unemp2_sorted <- unem2_gather[order(unem2_gather$Year),]
head(unemp2_sorted)
unem2 <- ts(unemp2_sorted$`unemployment rate`, start=c(2010,1), end=c(2019,12), frequency = 12)

# 2.1. DESCRIPTIVE ANALYSIS
descriptive_analysis(unem2, "Unemployment Rate (%)", "Year", "the Monthly U.S. Unemployment Rate for 2010-2019")
Acf(unem2, main = "ACF of the Monthly U.S. Unemployment Rate from 2010 to 2019")
mean(unem2)
sd(unem2)

# 2.2. LINEAR TREND MODEL
deterministic_model(unem2, model="Linear", title="the Monthly U.S. Unemployment Rate (2010-2019)", ylabel="Unemployment Rate (%)", xlabel="Year")

# 2.3. QUADRATIC TREND MODEL
deterministic_model(unem2, model="Quadratic", title="the Monthly U.S. Unemployment Rate (2010-2019)", ylabel="Unemployment Rate (%)", xlabel="Year")

# 2.4. STOCHASTIC TREND MODEL

# 2.4.1 Test of Normality and Stationarity
normality_test(unem2, "the U.S. Unemployment Rate") 
correlograms(unem2, "the U.S. Unemployment Rate")
adf.test(unem2)
unem2.transform = BoxCox.ar(unem2) # Data Transformation
unem2.transform$ci # lambda from 0.1 to 0.8
lambda = 0.45
unem2.BC = (unem2^lambda-1)/lambda
normality_test(unem2.BC, title = "the Time Series after Box-Cox Transformation")
plot(unem2.BC, ylab='Unemployment Rate (%)',xlab='Year',type='o', main = "Time Series Plot of the Monthly U.S. Unemployment Rate after Box-Cox Transformation")

# 2.4.2. Model Specification
differencing(unem2.BC)
diff2.unem2.BC = diff(unem2.BC,differences = 2) # Set data at Second differencing
plot(diff2.unem2.BC ,ylab='Unemployment Rate (%)',xlab='Year',type='o', main = "Time series plot of the Second Differenced Series (2010-2019)")
abline(h=0, col = 'red', lty = 2)

correlograms(diff2.unem2.BC, "the Second Differenced Series (2010-2019)")
eacf(diff2.unem2.BC,ar.max = 10, ma.max= 10)
bic_ms2 = armasubsets(y=diff2.unem2.BC,nar=5,nma=5,y.name='test',ar.method='ols')
plot(bic_ms2)

# 2.4.3. Model Fitting and Parameter Estimation
model2.021 <-  arima(unem2.BC, order=c(0,2,1), method='ML')
coeftest(model2.021)
model2.022 <-  arima(unem2.BC, order=c(0,2,2), method='ML')
coeftest(model2.022)
model2.122 <- arima(unem2.BC, order=c(1,2,2), method='ML')
coeftest(model2.122)
model2.323 <- arima(unem2.BC, order=c(3,2,3), method='ML')
coeftest(model2.323)
model2.523 <- arima(unem2.BC, order=c(5,2,3), method='ML')
coeftest(model2.523)
model2.525 <- arima(unem2.BC, order=c(5,2,5), method='ML')
coeftest(model2.525)

# 2.4.4. Residual Analysis
residual.analysis(model2.021)
residual.analysis(model2.022)
residual.analysis(model2.122)
residual.analysis(model2.323)
residual.analysis(model2.523)
residual.analysis(model2.525)

# 2.4.5. Model Selection (of the remaining candidate models)
sort.score(AIC(model2.122, model2.323, model2.523, model2.525), score = "aic")
sort.score(BIC(model2.122, model2.323, model2.523, model2.525), score = "bic")

AICc(model2.122)
AICc(model2.523)
AICc(model2.323)
AICc(model2.525)

# 2.4.6. Overfitting (of the best model)
model2.222 = arima(unem2.BC, order=c(2,2,2), method='ML')
coeftest(model2.222)
model2.123 = arima(unem2.BC, order=c(1,2,3), method='ML')
coeftest(model2.123)

# 2.5. FORECAST
fit_122 = Arima(unem2,c(1,2,2), lambda = 0.45)
predict_122 = forecast(fit_122, h=10)
predict_122
plot(predict_122, type = 'o', colour = "blue", xlab = "Year", ylab = "Unemployment Rate (%)", main = "Time Series Plot of the U.S. Unemployment Rate with a Ten-month Forecast")

# DISCUSSION
par(mfrow=c(1,2))
plot(predict_115,  type = 'o', colour = "blue", xlab = "Year", ylab = "Unemployment Rate (%)", main = "Ten-month Forecast (2000-2019)") 
plot(predict_122, type = 'o', colour = "blue", xlab = "Year", ylab = "Unemployment Rate (%)", main = "Ten-month Forecast (2010-2019)")
par(mfrow=c(1,1))
