miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)

jack_w_2023 <- sqrt(jack_u10_2023^2 + jack_v10_2023^2)
orla_w_2023 <- sqrt(orla_u10_2023^2 + orla_v10_2023^2)
fort_w_2023 <- sqrt(fort_u10_2023^2 + fort_v10_2023^2)

set.seed(5)

# generate a formula to use for finding coefficients
create_formula <- function(key_freqs) {
  sentence <- ""
  if (length(key_freqs)%%2 == 0){
    k <- length(key_freqs)/2 + 1
  } else {
    k <- round(length(key_freqs)/2, 0)
  }
  
  for (i in k:length(key_freqs)) {
    term <- paste("sin(2*pi*", key_freqs[i], "*t) + cos(2*pi*", key_freqs[i], "*t)")
    sentence <- paste(sentence, term, "+", sep = "")
  }
  return(as.formula(paste("z ~", substr(sentence, 1, nchar(sentence) - 1))))
}

t <- 1:length(miami_w_2023)
freq <- seq(from = -0.5, to = 0.5, length = length(miami_w_2023))

# fit a personal stl model to miami wind speed
# aim for an additive model

# fit the trend
# acf(miami_w_2023)

# quadratic
quadfit <- lm(miami_w_2023 ~ poly(t, 2, raw=TRUE))
z <- miami_w_2023 - quadfit$fitted.values
key_freqs <- c()
for(i in 1:length(my_periodogram(z))){
  if(my_periodogram(z)[i] > 170){
    key_freqs <- c(key_freqs, freq[i])
  }
}
key_freqs <- key_freqs[11:20]
seasonqfit <- lm(create_formula(key_freqs))
y <- z - seasonqfit$fitted.values
# fit arima model to y
y_model <- auto.arima(y, ic = "aicc")
# fit arma garch model to y 
yspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                    mean.model=list(armaOrder = c(1, 4), include.mean=TRUE))
y_garch <- ugarchfit(spec = yspec, data = y, solver="lbfgs")
y_garch@fit$solver$sol$par
y_gar_res <- y_garch@fit$residuals
y_armac <- y_garch@fit$coef[2:6]
y_garchc <- y_garch@fit$coef[7:9]

# plots
plot(miami_w_2023, type="l", main="Quadratic Trend")
lines(quadfit$fitted.values, col="red")

plot(freq, my_periodogram(z), type="l")

plot(z, type="l", main="Quadratic Trend with Seasonal Component")
lines(seasonqfit$fitted.values, col="red")

par(mfrow = c(2,2))
plot(y, type="l", main="Remainder with Quadratic Trend")
hist(y, breaks=100, freq=FALSE, main="Histogram")
lines(density(y), col="red")
qqnorm(y, main="QQ Plot")
acf(y)



# tampa wind speed
# quadratic
quadfit2 <- lm(tampa_w_2023 ~ poly(t, 2, raw=TRUE))
z2 <- tampa_w_2023 - quadfit2$fitted.values
key_freqs2 <- c()
for(i in 1:length(my_periodogram(z2))){
  if(my_periodogram(z2)[i] > 125){
    key_freqs2 <- c(key_freqs2, freq[i])
  }
}
key_freqs2 <- key_freqs2[11:20]
seasonqfit2 <- lm(create_formula(key_freqs2))
y2 <- z2 - seasonqfit2$fitted.values
# fit arima model to y2
y2_model <- auto.arima(y2, ic = "aicc")

y2spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                    mean.model=list(armaOrder = c(5, 2), include.mean=TRUE))
y2_garch <- ugarchfit(spec = y2spec, data = y2, solver="hybrid")
y2_garch@fit$solver$sol$par
y2_armac <- y2_garch@fit$solver$sol$par[2:8]
y2_garchc <- y2_garch@fit$solver$sol$par[9:11]
y2_gar_res <- y2_garch@fit$residuals

# plots
plot(tampa_w_2023, type="l", main="Quadratic Trend")
lines(quadfit2$fitted.values, col="red")

plot(freq, my_periodogram(z2), type="l")

plot(z2, type="l", main="Quadratic Trend with Seasonal Component")
lines(seasonqfit2$fitted.values, col="red")

par(mfrow = c(2,2))
plot(y2, type="l", main="Remainder with Quadratic Trend")
hist(y2, breaks=100, freq=FALSE, main="Histogram")
lines(density(y2), col="red")
qqnorm(y2, main="QQ Plot")
acf(y2)

# ARIMA RESIDUALS
# fit mvn distribution to residuals
res_mvn <- mvn("XXX", cbind(y_model$residuals, y2_model$residuals))

# fit ghyp distribution to residuals
res_ghyp <- stepAIC.ghyp(cbind(y_model$residuals, y2_model$residuals), silent = TRUE)

res1_ghyp <- stepAIC.ghyp(y_model$residuals, silent=TRUE)
res2_ghyp <- stepAIC.ghyp(y2_model$residuals, silent=TRUE)

# ljung box test
checkresiduals(y_model)
checkresiduals(y2_model)


# ARIMA forecasting
t <- 5832
y11 <- y[t] + y_model$coef[1]*y_model$residuals[t] + y_model$coef[2]*y_model$residuals[t-1] + y_model$coef[3]*y_model$residuals[t-2] + y_model$coef[4]*y_model$residuals[t-3]
y12 <- y11 + y_model$coef[2]*y_model$residuals[t] + y_model$coef[3]*y_model$residuals[t-1] + y_model$coef[4]*y_model$residuals[t-2]
y13 <- y12 + y_model$coef[3]*y_model$residuals[t] + y_model$coef[4]*y_model$residuals[t-1]
y14 <- y13 + y_model$coef[4]*y_model$residuals[t] 
y15 <- y14
y21 <- y2_model$coef[1]*y2[t] + y2_model$coef[2]*y2[t-1] + y2_model$coef[3]*y2[t-2] + y2_model$coef[4]*y2[t-3] + y2_model$coef[5]*y2[t-4] + y2_model$coef[6]*y2_model$residuals[t] + y2_model$coef[7]*y2_model$residuals[t-1]
y22 <- y2_model$coef[1]*y21 + y2_model$coef[2]*y2[t] + y2_model$coef[3]*y2[t-1] + y2_model$coef[4]*y2[t-2] + y2_model$coef[5]*y2[t-3] + y2_model$coef[7]*y2_model$residuals[t]
y23 <- y2_model$coef[1]*y22 + y2_model$coef[2]*y21 + y2_model$coef[3]*y2[t] + y2_model$coef[4]*y2[t-1] + y2_model$coef[5]*y2[t-2]
y24 <- y2_model$coef[1]*y23 + y2_model$coef[2]*y22 + y2_model$coef[3]*y21 + y2_model$coef[4]*y2[t] + y2_model$coef[5]*y2[t-1]
y25 <- y2_model$coef[1]*y24 + y2_model$coef[2]*y23 + y2_model$coef[3]*y22 + y2_model$coef[4]*y21 + y2_model$coef[5]*y2[t]

# ARIMA check forecasting inside time series lines up with what was observed
# consider after the first 1000 observations
k <- 1000
y11 <- y[k] + y_model$coef[1]*y_model$residuals[k] + y_model$coef[2]*y_model$residuals[k-1] + y_model$coef[3]*y_model$residuals[k-2] + y_model$coef[4]*y_model$residuals[k-3]
y12 <- y11 + y_model$coef[2]*y_model$residuals[k] + y_model$coef[3]*y_model$residuals[k-1] + y_model$coef[4]*y_model$residuals[k-2]
y13 <- y12 + y_model$coef[3]*y_model$residuals[k] + y_model$coef[4]*y_model$residuals[k-1]
y14 <- y13 + y_model$coef[4]*y_model$residuals[k] 
y15 <- y14
y21 <- y2_model$coef[1]*y2[k] + y2_model$coef[2]*y2[k-1] + y2_model$coef[3]*y2[k-2] + y2_model$coef[4]*y2[k-3] + y2_model$coef[5]*y2[k-4] + y2_model$coef[6]*y2_model$residuals[k] + y2_model$coef[7]*y2_model$residuals[k-1]
y22 <- y2_model$coef[1]*y21 + y2_model$coef[2]*y2[k] + y2_model$coef[3]*y2[k-1] + y2_model$coef[4]*y2[k-2] + y2_model$coef[5]*y2[k-3] + y2_model$coef[7]*y2_model$residuals[k]
y23 <- y2_model$coef[1]*y22 + y2_model$coef[2]*y21 + y2_model$coef[3]*y2[k] + y2_model$coef[4]*y2[k-1] + y2_model$coef[5]*y2[k-2]
y24 <- y2_model$coef[1]*y23 + y2_model$coef[2]*y22 + y2_model$coef[3]*y21 + y2_model$coef[4]*y2[k] + y2_model$coef[5]*y2[k-1]
y25 <- y2_model$coef[1]*y24 + y2_model$coef[2]*y23 + y2_model$coef[3]*y22 + y2_model$coef[4]*y21 + y2_model$coef[5]*y2[k]

y1preds <- rep(y15, times=495)
y1preds <- c(y11, y12, y13, y14, y15, y1preds)

y2preds <- c(y21, y22, y23, y24, y25, rep(0, times=495))
for(i in 6:500){
  y2preds[i] <- y2_model$coef[1]*y2preds[i-1] + y2_model$coef[2]*y2preds[i-2] + y2_model$coef[3]*y2preds[i-3] + y2_model$coef[4]*y2preds[i-4] + y2_model$coef[5]*y2preds[i-5]
}

yts <- 1:1500
plot(yts, y[1:1500], type="l", xlab="Index", ylab="Speed", main="Miami")
lines(yts[1001:1500], y1preds, col="blue")
abline(v = 1000, col="red", lty=2)

plot(yts, y2[1:1500], type="l", xlab="Index", ylab="Speed", main="Tampa")
lines(yts[1001:1500], y2preds, col="blue")
abline(v = 1000, col="red", lty=2)


# ARMA GARCH

# ARMA GARCH residuals distribution
garch_res_ghyp <- stepAIC.ghyp(cbind(y_garch@fit$residuals, y2_garch@fit$residuals), silent = TRUE)
garch_res_ghyp$best.model

# ARMA GARCH forecasting
t <- 5832

y11g <- y_armac[1]*y[t] + y_armac[2]*y_gar_res[t] + y_armac[3]*y_gar_res[t-1] + y_armac[4]*y_gar_res[t-2] + y_armac[5]*y_gar_res[t-3]
y12g <- y_armac[1]*y11 + y_armac[3]*y_gar_res[t] + y_armac[4]*y_gar_res[t-1] + y_armac[5]*y_gar_res[t-2]
y13g <- y_armac[1]*y12 + y_armac[4]*y_gar_res[t] + y_armac[5]*y_gar_res[t-1]
y14g <- y_armac[1]*y13 + y_armac[5]*y_gar_res[t] 
y15g <- y_armac[1]*y14
y21g <- y2_armac[1]*y2[t] + y2_armac[2]*y2[t-1] + y2_armac[3]*y2[t-2] + y2_armac[4]*y2[t-3] + y2_armac[5]*y2[t-4] + y2_armac[6]*y2_gar_res[t] + y2_armac[7]*y2_gar_res[t-1]
y22g <- y2_armac[1]*y21 + y2_armac[2]*y2[t] + y2_armac[3]*y2[t-1] + y2_armac[4]*y2[t-2] + y2_armac[5]*y2[t-3] + y2_armac[7]*y2_gar_res[t]
y23g <- y2_armac[1]*y22 + y2_armac[2]*y21 + y2_armac[3]*y2[t] + y2_armac[4]*y2[t-1] + y2_armac[5]*y2[t-2]
y24g <- y2_armac[1]*y23 + y2_armac[2]*y22 + y2_armac[3]*y21 + y2_armac[4]*y2[t] + y2_armac[5]*y2[t-1]
y25g <- y2_armac[1]*y24 + y2_armac[2]*y23 + y2_armac[3]*y22 + y2_armac[4]*y21 + y2_armac[5]*y2[t]

y1predsg <- c(y11g, y12g, y13g, y14g, y15g, rep(0, times=495))
for(i in 6:500){
  y1predsg[i] <- y_armac[1]*y1predsg[i-1]
}

y2predsg <- c(y21g, y22g, y23g, y24g, y25g, rep(0, times=495))
for(i in 6:500){
  y2predsg[i] <- y2_armac[1]*y2preds[i-1] + y2_armac[2]*y2preds[i-2] + y2_armac[3]*y2preds[i-3] + y2_armac[4]*y2preds[i-4] + y2_armac[5]*y2preds[i-5]
}

yts <- 1:1500
plot(yts, y[1:1500], type="l", xlab="Index", ylab="Speed", main="Miami")
lines(yts[1001:1500], y1predsg, col="blue")
abline(v = 1000, col="red", lty=2)

plot(yts, y2[1:1500], type="l", xlab="Index", ylab="Speed", main="Tampa")
lines(yts[1001:1500], y2preds, col="blue")
abline(v = 1000, col="red", lty=2)



#ARMA GARCH conditional variance forecasting
sigma1 <- coredata(sigma(y_garch))[t,1]
sigma2 <- coredata(sigma(y2_garch))[t,1]

s11 <- y_garchc[1] + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
s12 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3])) + (y_garchc[2] + y_garchc[3])*(y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
s13 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^2
s14 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2 + (y_garchc[2] + y_garchc[3])^3) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^3
s15 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2 + (y_garchc[2] + y_garchc[3])^3 + (y_garchc[2] + y_garchc[3])^4) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^4
s21 <- y2_garchc[1] + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
s22 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3])) + (y2_garchc[2] + y2_garchc[3])*(y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
s23 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^2
s24 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2 + (y2_garchc[2] + y2_garchc[3])^3) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^3
s25 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2 + (y2_garchc[2] + y2_garchc[3])^3 + (y2_garchc[2] + y2_garchc[3])^4) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^4

# check forecasting inside time series lines up with what was observed
k <- 1000

y11g <- y_armac[1]*y[k] + y_armac[2]*y_gar_res[k] + y_armac[3]*y_gar_res[k-1] + y_armac[4]*y_gar_res[k-2] + y_armac[5]*y_gar_res[k-3]
y12g <- y_armac[1]*y11 + y_armac[3]*y_gar_res[k] + y_armac[4]*y_gar_res[k-1] + y_armac[5]*y_gar_res[k-2]
y13g <- y_armac[1]*y12 + y_armac[4]*y_gar_res[k] + y_armac[5]*y_gar_res[k-1]
y14g <- y_armac[1]*y13 + y_armac[5]*y_gar_res[k] 
y15g <- y_armac[1]*y14
y21g <- y2_armac[1]*y2[k] + y2_armac[2]*y2[k-1] + y2_armac[3]*y2[k-2] + y2_armac[4]*y2[k-3] + y2_armac[5]*y2[k-4] + y2_armac[6]*y2_gar_res[k] + y2_armac[7]*y2_gar_res[k-1]
y22g <- y2_armac[1]*y21 + y2_armac[2]*y2[k] + y2_armac[3]*y2[k-1] + y2_armac[4]*y2[k-2] + y2_armac[5]*y2[k-3] + y2_armac[7]*y2_gar_res[k]
y23g <- y2_armac[1]*y22 + y2_armac[2]*y21 + y2_armac[3]*y2[k] + y2_armac[4]*y2[k-1] + y2_armac[5]*y2[k-2]
y24g <- y2_armac[1]*y23 + y2_armac[2]*y22 + y2_armac[3]*y21 + y2_armac[4]*y2[k] + y2_armac[5]*y2[k-1]
y25g <- y2_armac[1]*y24 + y2_armac[2]*y23 + y2_armac[3]*y22 + y2_armac[4]*y21 + y2_armac[5]*y2[k]


# conditional variance formula
miami_s_fc <- function(tau){
  t <- 5832
  summand <- rep(0, times = tau)
  summand[1] <- 1
  if(tau > 1){
    for(i in 2:tau){
      summand[i] <- (y_garchc[2] + y_garchc[3])^(i-1)
    }
  }
  sig <- y_garchc[1]*sum(summand) + (y_garchc[2] + y_garchc[3])^(tau-1)*(y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
  return(sig)
}

tampa_s_fc <- function(tau){
  t <- 5832
  summand <- rep(0, times = tau)
  summand[1] <- 1
  if(tau > 1){
    for(i in 2:tau){
      summand[i] <- (y2_garchc[2] + y2_garchc[3])^(i-1)
    }
  }
  sig <- y2_garchc[1]*sum(summand) + (y2_garchc[2] + y2_garchc[3])^(tau-1)*(y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
  return(sig)
}

xtss <- 1:(length(miami_w_2023)+48)
miami_s <- rep(0, times = 48)
tampa_s <- rep(0, times = 48)
for(i in 1:48){
  miami_s[i] <- miami_s_fc(i)
  tampa_s[i] <- tampa_s_fc(i)
}

# plot showing last 1 week and then prediction for next week
plot(xtss[(t-168):length(xtss)], c(coredata(sigma(y_garch))[(t-168):t,1], miami_s), 
     type="l", xlab="Index", ylab="Variance", main="Miami")
abline(v = length(miami_w_2023), col="red", lty=2)

plot(xtss[(t-168):length(xtss)], c(coredata(sigma(y2_garch))[(t-168):t,1], tampa_s), 
     type="l", xlab="Index", ylab="Variance", main="Tampa")
abline(v = length(tampa_w_2023), col="red", lty=2)





# VAR model
ys_var <- VAR(cbind(y, y2), p=1, type="none")
# VAR prediction
phi_var <- matrix(c(ys_var$varresult$y$coefficients, ys_var$varresult$y2$coefficients),
                     nrow=2, ncol=2, byrow=TRUE)
Y_fin <- c(y[5832], y2[5832])

var_forecast <- function(tau){
  fore <- Y_fin
  for (i in 1:tau){
    fore <- phi_var %*% fore
  }
  return(fore)
}

# create a plot showing the time series and then the forecasts
miamiforecast <- c(rep(0, times = 504))
tampaforecast <- c(rep(0, times = 504))

xts <- 1:(length(miami_w_2023)+504)

for(i in 1:504){
  miamiforecast[i] <- var_forecast(i)[1]
  tampaforecast[i] <- var_forecast(i)[2]
}

par(mfrow=c(1,2))
plot(xts[5000:length(xts)], c(y, miamiforecast)[5000:length(xts)], 
     type="l", main="Miami", xlab="Index", ylab="Speed")
abline(v = length(miami_w_2023), col="red", lty=2)

plot(xts[5000:length(xts)], c(y2, tampaforecast)[5000:length(xts)], 
     type="l", main="Tampa", xlab="Index", ylab="Speed")
abline(v = length(tampa_w_2023), col="red", lty=2)

y_var_res <- ys_var$varresult$y$residuals
y2_var_res <- ys_var$varresult$y2$residuals

var_res_mvn <- mvn("XXX", cbind(y_var_res, y2_var_res))
var_res_ghyp <- stepAIC.ghyp(cbind(y_var_res, y2_var_res), silent=TRUE)



# tallahassee ARMA-GARCH model
quadfit3 <- lm(talla_w_2023 ~ poly(t, 2, raw=TRUE))
z3 <- talla_w_2023 - quadfit3$fitted.values
key_freqs3 <- c()
for(i in 1:length(my_periodogram(z3))){
  if(my_periodogram(z3)[i] > 31.5){
    key_freqs3 <- c(key_freqs3, freq[i])
  }
}
key_freqs3 <- key_freqs3[11:20]
seasonqfit3 <- lm(create_formula(key_freqs3))
y3 <- z3 - seasonqfit3$fitted.values
# fit arima model to y3
y3_model <- auto.arima(y3, ic = "aicc")
# fit arma garch model to y3 using arima choice 
y3spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(1,1), include.mean=TRUE))
y3_garch <- ugarchfit(spec = y3spec, data = y3, solver="hybrid")
y3_garch@fit$solver$sol$par
y3_armac <- y3_garch@fit$solver$sol$par[2:3]
y3_garchc <- y3_garch@fit$solver$sol$par[4:6]
y3_gar_res <- y3_garch@fit$residuals

# miami, talla
garch2_res_ghyp <- stepAIC.ghyp(cbind(y_garch@fit$residuals, y3_garch@fit$residuals), silent = TRUE)
garch2_res_ghyp$best.model

# tampa, talla
garch3_res_ghyp <- stepAIC.ghyp(cbind(y2_garch@fit$residuals, y3_garch@fit$residuals), silent = TRUE)
garch3_res_ghyp$best.model

# jackson
quadfit4 <- lm(jack_w_2023 ~ poly(t, 2, raw=TRUE))
z4 <- jack_w_2023 - quadfit4$fitted.values
key_freqs4 <- c()
for(i in 1:length(my_periodogram(z4))){
  if(my_periodogram(z4)[i] > 47){
    key_freqs4 <- c(key_freqs4, freq[i])
  }
}
key_freqs4 <- key_freqs4[11:20]
seasonqfit4 <- lm(create_formula(key_freqs4))
y4 <- z4 - seasonqfit4$fitted.values
# fit arima model to y4
y4_model <- auto.arima(y4, ic = "aicc")
# fit arma garch model to y4 using arima choice 
y4spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(2,1), include.mean=TRUE))
y4_garch <- ugarchfit(spec = y4spec, data = y4, solver="hybrid")
y4_garch@fit$solver$sol$par
y4_gar_res <- y4_garch@fit$residuals

# orlando
quadfit5 <- lm(orla_w_2023 ~ poly(t, 2, raw=TRUE))
z5 <- orla_w_2023 - quadfit5$fitted.values
key_freqs5 <- c()
for(i in 1:length(my_periodogram(z5))){
  if(my_periodogram(z5)[i] > 77){
    key_freqs5 <- c(key_freqs5, freq[i])
  }
}
key_freqs5 <- key_freqs5[11:20]
seasonqfit5 <- lm(create_formula(key_freqs5))
y5 <- z5 - seasonqfit5$fitted.values
# fit arima model to y5
y5_model <- auto.arima(y5, ic = "aicc")
# fit arma garch model to y5 using arima choice 
y5spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(1,1), include.mean=TRUE))
y5_garch <- ugarchfit(spec = y5spec, data = y5, solver="hybrid")
y5_garch@fit$solver$sol$par
y5_gar_res <- y5_garch@fit$residuals

# fort meyers
quadfit6 <- lm(fort_w_2023 ~ poly(t, 2, raw=TRUE))
z6 <- fort_w_2023 - quadfit6$fitted.values
key_freqs6 <- c()
for(i in 1:length(my_periodogram(z6))){
  if(my_periodogram(z6)[i] > 97.75){
    key_freqs6 <- c(key_freqs6, freq[i])
  }
}
key_freqs6 <- key_freqs6[11:20]
seasonqfit6 <- lm(create_formula(key_freqs6))
y6 <- z6 - seasonqfit6$fitted.values
# fit arima model to y6
y6_model <- auto.arima(y6, ic = "aicc")
# fit arma garch model to y6 using arima choice 
y6spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(3,1), include.mean=TRUE))
y6_garch <- ugarchfit(spec = y6spec, data = y6, solver="hybrid")
y6_garch@fit$solver$sol$par
y6_gar_res <- y6_garch@fit$residuals




# try fitting another seasonal component to the residuals
key_freqs11 <- c()
for(i in 1:length(my_periodogram(y_gar_res))){
  if(my_periodogram(y_gar_res)[i] > 1.9){
    key_freqs11 <- c(key_freqs11, freq[i])
  }
}
key_freqs11 <- key_freqs11[12:21]
seasongar1fit <- lm(create_formula(key_freqs11))
y_gar_res_new <- y_gar_res - seasongar1fit$fitted.values

key_freqs22 <- c()
for(i in 1:length(my_periodogram(y2_gar_res))){
  if(my_periodogram(y2_gar_res)[i] > 1.8){
    key_freqs22 <- c(key_freqs22, freq[i])
  }
}
key_freqs22 <- key_freqs22[11:20]
seasongar2fit <- lm(create_formula(key_freqs22))
y2_gar_res_new <- y2_gar_res - seasongar2fit$fitted.values
