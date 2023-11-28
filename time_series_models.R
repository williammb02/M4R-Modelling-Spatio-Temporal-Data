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
plot(z, type="l")
