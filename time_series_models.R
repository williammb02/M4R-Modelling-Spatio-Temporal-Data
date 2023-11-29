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

# cubic trend
polyfit <- lm(miami_w_2023 ~ poly(t, 3, raw=TRUE))
plot(miami_w_2023, type="l", main="Cubic Trend")
lines(polyfit$fitted.values, col="red")

z <- miami_w_2023 - polyfit$fitted.values
plot(freq, my_periodogram(z), type="l")

# identify key periods in the data 
key_freqs <- c()
powers <- c()
for(i in 1:length(my_periodogram(z))){
  if(my_periodogram(z)[i] > 135){
    key_freqs <- c(key_freqs, freq[i])
    powers <- c(powers, my_periodogram(z)[i])
  }
}

seasonfit_test <- lm(create_formula(key_freqs))
plot(z, type="l", main="Cubic Trend with Seasonal Component of 28 Key Frequencies")
lines(seasonfit_test$fitted.values, col="red")

y <- z - seasonfit_test$fitted.values
par(mfrow = c(2,2))
plot(y, type="l", main="Remainder with Cubic Trend, 28 Frequencies")
hist(y, breaks=100, freq=FALSE)
lines(density(y), col="red")
qqnorm(y)
acf(y)

# try with new threshold
key_freqs2 <- c()
powers2 <- c()
for(i in 1:length(my_periodogram(z))){
  if(my_periodogram(z)[i] > 10){
    key_freqs2 <- c(key_freqs2, freq[i])
    powers2 <- c(powers2, my_periodogram(z)[i])
  }
}

seasonfit <- lm(create_formula(key_freqs2))
plot(z, type="l", main="Cubic Trend with Seasonal Component of 296 Key Frequencies")
lines(seasonfit$fitted.values, col="red")

y2 <- z - seasonfit$fitted.values
par(mfrow = c(2,2))
plot(y2, type="l", main="Remainder with Cubic Trend, 296 Frequencies")
hist(y2, breaks=100, freq=FALSE)
lines(density(y2), col="red")
qqnorm(y2)
acf(y2)

# quadratic
# need to run all code here separate to the cubic as we redefine z
quadfit <- lm(miami_w_2023 ~ poly(t, 2, raw=TRUE))
plot(miami_w_2023, type="l", main="Quadratic Trend")
lines(quadfit$fitted.values, col="red")
z <- miami_w_2023 - quadfit$fitted.values
plot(freq, my_periodogram(z), type="l")
key_freqs3 <- c()
powers3 <- c()
for(i in 1:length(my_periodogram(z))){
  if(my_periodogram(z)[i] > 10){
    key_freqs3 <- c(key_freqs3, freq[i])
    powers3 <- c(powers3, my_periodogram(z)[i])
  }
}
seasonqfit <- lm(create_formula(key_freqs3))
plot(z, type="l", main="Quadratic Trend with Seasonal Component of 298 Key Frequencies")
lines(seasonqfit$fitted.values, col="red")
y3 <- z - seasonqfit$fitted.values
par(mfrow = c(2,2))
plot(y3, type="l", main="Remainder with Quadratic Trend, 298 Frequencies")
hist(y3, breaks=100, freq=FALSE, main="Histogram")
lines(density(y3), col="red")
qqnorm(y3, main="QQ Plot")
acf(y3)

# smooth periodogram
# does not create a good fit
spz <- smooth.periodogram(z)
spz_vals <- spz$smooth.periodogram
spz_freq <- spz$lambda
key_freqs4 <- c()
for(i in 1:length(spz_vals)){
  if(spz_vals[i] > 2){
    key_freqs4 <- c(key_freqs4, spz_freq[i])
  }
}
seasonqfit2 <- lm(create_formula(key_freqs4))
plot(z, type="l", main="Quadratic Trend with Seasonal Component of 129 Key Frequencies")
lines(seasonqfit2$fitted.values, col="red")
y4 <- z - seasonqfit2$fitted.values
par(mfrow = c(2,2))
plot(y4, type="l", main="Remainder with Quadratic Trend, 129 Frequencies")
hist(y4, breaks=100, freq=FALSE, main="Histogram")
lines(density(y4), col="red")
qqnorm(y4, main="QQ Plot")
acf(y4)


# try differencing?
freq <- seq(from = -0.5, to = 0.5, length = length(diff(miami_w_2023)))

plot(diff(miami_w_2023), type="l")
pacf(diff(miami_w_2023), lag.max = 1000)
plot(freq, my_periodogram(diff(miami_w_2023)), type="l")

key_freqs <- c()
powers <- c()
for(i in 1:length(my_periodogram(diff(miami_w_2023)))){
  if(my_periodogram(diff(miami_w_2023))[i] > 2){
    key_freqs <- c(key_freqs, freq[i])
    powers <- c(powers, my_periodogram(diff(miami_w_2023))[i])
  }
}
key_periods = 1/key_freqs
round(key_periods, 0)
