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
acf(miami_w_2023)

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
                    mean.model=list(armaOrder = c(1, 1), include.mean=TRUE))
y_garch <- ugarchfit(spec = yspec, data = y)
y_garch@fit$solver$sol$pars

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
                    mean.model=list(armaOrder = c(1, 1), include.mean=TRUE))
y2_garch <- ugarchfit(spec = y2spec, data = y2)
y2_garch@fit$solver$sol$pars

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


# fit mvn distribution to residuals
res_mvn <- mvn("XXX", cbind(y_model$residuals, y2_model$residuals))

# fit ghyp distribution to residuals
res_ghyp <- stepAIC.ghyp(cbind(y_model$residuals, y2_model$residuals), silent = TRUE)

res1_ghyp <- stepAIC.ghyp(y_model$residuals, silent=TRUE)
res2_ghyp <- stepAIC.ghyp(y2_model$residuals, silent=TRUE)

# ljung box test
checkresiduals(y_model)
checkresiduals(y2_model)



# VAR model
ys_var <- VAR(cbind(y, y2), p=1, type="none")

y_var_res <- ys_var$varresult$y$residuals
y2_var_res <- ys_var$varresult$y2$residuals


var_res_mvn <- mvn("XXX", cbind(y_var_res, y2_var_res))
var_res_ghyp <- stepAIC.ghyp(cbind(y_var_res, y2_var_res), silent=TRUE)
