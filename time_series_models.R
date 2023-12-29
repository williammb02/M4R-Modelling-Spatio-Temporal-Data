miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)

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


#ARIMA FORECASTING
t <- 5832
y11 <- miami_w_2023[t] + y_model$coef[1]*y_model$residuals[t] + y_model$coef[2]*y_model$residuals[t-1] + y_model$coef[3]*y_model$residuals[t-2] + y_model$coef[4]*y_model$residuals[t-3]
y12 <- y11 + y_model$coef[2]*y_model$residuals[t] + y_model$coef[3]*y_model$residuals[t-1] + y_model$coef[4]*y_model$residuals[t-2]
y13 <- y12 + y_model$coef[3]*y_model$residuals[t] + y_model$coef[4]*y_model$residuals[t-1]
y14 <- y13 + y_model$coef[4]*y_model$residuals[t] 
y15 <- y14
y21 <- y2_model$coef[1]*tampa_w_2023[t] + y2_model$coef[2]*tampa_w_2023[t-1] + y2_model$coef[3]*tampa_w_2023[t-2] + y2_model$coef[4]*tampa_w_2023[t-3] + y2_model$coef[5]*tampa_w_2023[t-4] + y2_model$coef[6]*y2_model$residuals[t] + y2_model$coef[7]*y2_model$residuals[t-1]
y22 <- y2_model$coef[1]*y21 + y2_model$coef[2]*tampa_w_2023[t] + y2_model$coef[3]*tampa_w_2023[t-1] + y2_model$coef[4]*tampa_w_2023[t-2] + y2_model$coef[5]*tampa_w_2023[t-3] + y2_model$coef[7]*y2_model$residuals[t]
y23 <- y2_model$coef[1]*y22 + y2_model$coef[2]*y21 + y2_model$coef[3]*tampa_w_2023[t] + y2_model$coef[4]*tampa_w_2023[t-1] + y2_model$coef[5]*tampa_w_2023[t-2]
y24 <- y2_model$coef[1]*y23 + y2_model$coef[2]*y22 + y2_model$coef[3]*y21 + y2_model$coef[4]*tampa_w_2023[t] + y2_model$coef[5]*tampa_w_2023[t-1]
y25 <- y2_model$coef[1]*y24 + y2_model$coef[2]*y23 + y2_model$coef[3]*y22 + y2_model$coef[4]*y21 + y2_model$coef[5]*tampa_w_2023[t]




# ARMA GARCH RESIDUALS
garch_res_ghyp <- stepAIC.ghyp(cbind(y_garch@fit$residuals, y2_garch@fit$residuals), silent = TRUE)
garch_res_ghyp$best.model

t <- 5832

y11g <- y_armac[1]*miami_w_2023[t] + y_armac[2]*y_gar_res[t] + y_armac[3]*y_gar_res[t-1] + y_armac[4]*y_gar_res[t-2] + y_armac[5]*y_gar_res[t-3]
y12g <- y_armac[1]*y11 + y_armac[3]*y_gar_res[t] + y_armac[4]*y_gar_res[t-1] + y_armac[5]*y_gar_res[t-2]
y13g <- y_armac[1]*y12 + y_armac[4]*y_gar_res[t] + y_armac[5]*y_gar_res[t-1]
y14g <- y_armac[1]*y13 + y_armac[5]*y_gar_res[t] 
y15g <- y_armac[1]*y14
y21g <- y2_armac[1]*tampa_w_2023[t] + y2_armac[2]*tampa_w_2023[t-1] + y2_armac[3]*tampa_w_2023[t-2] + y2_armac[4]*tampa_w_2023[t-3] + y2_armac[5]*tampa_w_2023[t-4] + y2_armac[6]*y2_gar_res[t] + y2_armac[7]*y2_gar_res[t-1]
y22g <- y2_armac[1]*y21 + y2_armac[2]*tampa_w_2023[t] + y2_armac[3]*tampa_w_2023[t-1] + y2_armac[4]*tampa_w_2023[t-2] + y2_armac[5]*tampa_w_2023[t-3] + y2_armac[7]*y2_gar_res[t]
y23g <- y2_armac[1]*y22 + y2_armac[2]*y21 + y2_armac[3]*tampa_w_2023[t] + y2_armac[4]*tampa_w_2023[t-1] + y2_armac[5]*tampa_w_2023[t-2]
y24g <- y2_armac[1]*y23 + y2_armac[2]*y22 + y2_armac[3]*y21 + y2_armac[4]*tampa_w_2023[t] + y2_armac[5]*tampa_w_2023[t-1]
y25g <- y2_armac[1]*y24 + y2_armac[2]*y23 + y2_armac[3]*y22 + y2_armac[4]*y21 + y2_armac[5]*tampa_w_2023[t]

#ARMA GARCH conditional variances
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




# VAR model
ys_var <- VAR(cbind(y, y2), p=1, type="none")
# prediction
phi_var <- matrix(c(ys_var$varresult$y$coefficients, ys_var$varresult$y2$coefficients),
                     nrow=2, ncol=2, byrow=TRUE)
Y_fin <- c(miami_w_2023[5832], tampa_w_2023[5832])

var_forecast <- function(tau){
  fore <- Y_fin
  for (i in 1:tau){
    fore <- phi_var %*% fore
  }
  return(fore)
}

y_var_res <- ys_var$varresult$y$residuals
y2_var_res <- ys_var$varresult$y2$residuals


var_res_mvn <- mvn("XXX", cbind(y_var_res, y2_var_res))
var_res_ghyp <- stepAIC.ghyp(cbind(y_var_res, y2_var_res), silent=TRUE)



# TAIL DEPENDENCE
tde_lag <- function(u, tau){
  q1 <- quantile(miami_w_2023, u)
  q2 <- quantile(tampa_w_2023, u)
  num <- 0
  den <- 0
  l <- length(miami_w_2023)
  for (i in 1:(l-tau)){
    if ((miami_w_2023[i] > u) & (tampa_w_2023[i+tau] > u)){
      num <- num + 1
    }
    if (miami_w_2023[i] > u){
      den <- den + 1
    }
  }
  return(c(num/l, den/l, num/den))
}

tde_lag(0.8, 0)
tde_lag(0.9, 0)
tde_lag(0.95, 0)

tde_lag(0.95, 12)
tde_lag(0.95, 24)
tde_lag(0.95, 168)

tde_lag_2 <- function(x, y, u, tau){
  q1 <- quantile(x, u)
  q2 <- quantile(y, u)
  num <- 0
  den <- 0
  l <- length(x)
  for (i in 1:(l-tau)){
    if ((x[i] > u) & (y[i+tau] > u)){
      num <- num + 1
    }
    if (x[i] > u){
      den <- den + 1
    }
  }
  return(c(num/l, den/l, num/den))
}

tde_lag_2(y_model$residuals, y2_model$residuals, 0.8, 0)
tde_lag_2(y_model$residuals, y2_model$residuals, 0.9, 0)
tde_lag_2(y_model$residuals, y2_model$residuals, 0.95, 0)

tde_lag_2(y_gar_res, y2_gar_res, 0.8, 0)
tde_lag_2(y_gar_res, y2_gar_res, 0.9, 0)
tde_lag_2(y_gar_res, y2_gar_res, 0.95, 0)

tde_lag_2(y_var_res, y2_var_res, 0.8, 0)
tde_lag_2(y_var_res, y2_var_res, 0.9, 0)
tde_lag_2(y_var_res, y2_var_res, 0.95, 0)
