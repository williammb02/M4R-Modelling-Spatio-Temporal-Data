miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)


# fit a personal stl model to miami wind speed
# aim for an additive model

# fit the trend
acf(miami_w_2023)

t <- 1:length(miami_w_2023)
polyfit <- lm(miami_w_2023 ~ poly(t, 3, raw=TRUE))

plot(miami_w_2023, type="l")
lines(polyfit$fitted.values, col="red")

z <- miami_w_2023 - polyfit$fitted.values
# plot(z, type="l")
freq <- seq(from = -0.5, to = 0.5, length = length(miami_w_2023))
plot(freq, periodogram(z), type="l")

# identify key periods in the data 
key_freqs <- c()
powers <- c()
for(i in 1:length(periodogram(z))){
  if(periodogram(z)[i] > 135){
    key_freqs <- c(key_freqs, freq[i])
    powers <- c(powers, periodogram(z)[i])
  }
}

key_periods = 1/key_freqs
round(key_periods[15:28], 0)

# consider using all of them - index 15:28
create_formula <- function(key_freqs) {
  sentence <- ""
  if (length(key_freqs)%%2 == 0){
    k <- length(key_freqs)/2 + 1
  } else {
    k <- round(length(key_freqs), 0)
  }
  
  for (i in k:length(key_freqs)) {
    term <- paste("sin(2*pi*", key_freqs[i], "*t) + cos(2*pi*", key_freqs[i], "*t)")
    sentence <- paste(sentence, term, "+", sep = "")
  }
  return(as.formula(paste("z ~", substr(sentence, 1, nchar(sentence) - 1))))
}
create_formula(key_freqs)


seasonfit_test <- lm(create_formula(key_freqs))

plot(z, type="l")
lines(seasonfit_test$fitted.values, col="blue")

y <- z - seasonfit_test$fitted.values
plot(y, type="l")

# try with new threshold
key_freqs2 <- c()
powers2 <- c()
for(i in 1:length(periodogram(z))){
  if(periodogram(z)[i] > 5){
    key_freqs2 <- c(key_freqs2, freq[i])
    powers2 <- c(powers2, periodogram(z)[i])
  }
}

seasonfit <- lm(create_formula(key_freqs2))
plot(z, type="l")
lines(seasonfit$fitted.values, col="red")

y2 <- z - seasonfit$fitted.values
plot(periodogram(y2), type="l")

# try differencing?
freq <- seq(from = -0.5, to = 0.5, length = length(diff(miami_w_2023)))

plot(diff(miami_w_2023), type="l")
pacf(diff(miami_w_2023), lag.max = 1000)
plot(freq, periodogram(diff(miami_w_2023)), type="l")

key_freqs <- c()
powers <- c()
for(i in 1:length(periodogram(diff(miami_w_2023)))){
  if(periodogram(diff(miami_w_2023))[i] > 2){
    key_freqs <- c(key_freqs, freq[i])
    powers <- c(powers, periodogram(diff(miami_w_2023))[i])
  }
}
key_periods = 1/key_freqs
round(key_periods, 0)
