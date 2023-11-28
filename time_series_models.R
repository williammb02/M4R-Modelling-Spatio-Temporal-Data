miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)


# fit a personal stl model to miami wind speed
# aim for an additive model

# fit the trend
acf(miami_w_2023)

polyfit <- lm(miami_w_2023 ~ poly(1:length(miami_w_2023), 3, raw=TRUE))
# plot(miami_w_2023, type="l")
# lines(polyfit$fitted.values)
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
# pick index 28, 27, 21, 17





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
