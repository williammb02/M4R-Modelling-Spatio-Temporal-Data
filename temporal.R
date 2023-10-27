# import precipitation and wind for miami, tampa, and tallahassee
miami_tp <- import_climate_data("tp", miami_index[1], miami_index[2])
tampa_tp <- import_climate_data("tp", tampa_index[1], tampa_index[2])
tallahassee_tp <- import_climate_data("tp", tallahassee_index[1], tallahassee_index[2])

miami_u10 <- import_climate_data("u10", miami_index[1], miami_index[2])
tampa_u10 <- import_climate_data("u10", tampa_index[1], tampa_index[2])
tallahassee_u10 <- import_climate_data("u10", tallahassee_index[1], tallahassee_index[2])

miami_v10 <- import_climate_data("v10", miami_index[1], miami_index[2])
tampa_v10 <- import_climate_data("v10", tampa_index[1], tampa_index[2])
tallahassee_v10 <- import_climate_data("v10", tallahassee_index[1], tallahassee_index[2])

# create scatter plots for wind and precipitation at 3 locations in Florida
par(mfrow = c(2,3))
plot(miami_u10, miami_tp, main = "Miami", xlab="u10 (m s^-1)", ylab="Precipitation (m)")
plot(tampa_u10, tampa_tp, main = "Tampa", xlab="u10 (m s^-1)", ylab="Precipitation (m)")
plot(tallahassee_u10, tallahassee_tp, main = "Tallahassee", xlab="u10 (m s^-1)", ylab="Precipitation (m)")
plot(miami_v10, miami_tp, main = "Miami", xlab="v10 (m s^-1)", ylab="Precipitation (m)")
plot(tampa_v10, tampa_tp, main = "Tampa", xlab="v10 (m s^-1)", ylab="Precipitation (m)")
plot(tallahassee_v10, tallahassee_tp, main = "Tallahassee", xlab="v10 (m s^-1)", ylab="Precipitation (m)")

# create periodograms for the data
par(mfrow = c(1,3))
freq <- seq(-0.5, 0.5, length.out = length(miami_tp))

plot(freq, periodogram(miami_tp), type = "l", xlab = "Frequency", ylab = "Periodogram", main = "Miami")
plot(freq, periodogram(tampa_tp), type = "l", xlab = "Frequency", ylab = "Periodogram", main = "Tampa")
plot(freq, periodogram(tallahassee_tp), type = "l", xlab = "Frequency", ylab = "Periodogram", main = "Tallahassee")

# look at first differences of Miami precipitation
par(mfrow = c(1,2))
plot(miami_tp, type = "l", main = "Miami Total Precipitation")
acf(miami_tp, main = "Empirical Autocorrelation")

# consider NYC total precipitation to see if trend differs
nyc_tp <- import_climate_data("tp", nyc_index[1], nyc_index[2])
ts_nyc_tp <- ts(nyc_tp, frequency = 12)
nyc_stl <- stl(ts_nyc_tp, s.window = "periodic", t.window = length(ts_nyc_tp))
plot(nyc_stl, main = "STL Decomposition of NYC Total Precipitation")

nyc_remainder <- nyc_stl$time.series[,3]
par(mfrow = c(1,3))
acf(nyc_remainder)
hist(nyc_remainder)
plot(periodogram(as.vector(nyc_remainder)), type = "l")

# stl decomposition of precipitation
ts_miami_tp <- ts(miami_tp, frequency = 12)
miami_tp_trends <- stl(ts_miami_tp, s.window = "periodic", t.window = length(ts_miami_tp))
plot(miami_tp_trends, main = "Additive STL Decomposition of Total Precipitation in Miami")

# exploratory plots for the remainder
miami_remainder <- miami_tp_trends$time.series[,3]
par(mfrow = c(2,2))
# change between acf and pacf for an appropriate plot, pacf removes irrelevant 1.0 term
pacf(miami_remainder)
hist(miami_remainder, main = "Histogram")
plot(periodogram(as.vector(miami_remainder)), type = "l", ylab = "", main = "Periodogram")
qqnorm(miami_remainder)
qqline(miami_remainder, col = "red")

# fit ARIMA model to Miami total precipitation
miami_model <- auto.arima(miami_tp, ic = "aicc")
print(miami_model)
checkresiduals(miami_model)

# fit an ARIMA model to the remainders
miami_remainder_model <- auto.arima(miami_remainder, ic = "aicc")
print(miami_remainder_model)
checkresiduals(miami_remainder_model)

# fit a generalised hyperbolic distribution to the data
miami_rem_ghyp <- stepAIC.ghyp(miami_remainder, silent = TRUE)
hist(miami_rem_ghyp$best.model, main = "Histogram")

# analysis of absolute value or square of residuals
par(mfrow = c(1,2))
pacf(abs(miami_remainder))
pacf(miami_remainder^2)

# try stl decomposition of total precipitation, under first differencing
ts_miami_tp <- ts(diff(miami_tp), frequency = 12)
miami_tp_trends <- stl(ts_miami_tp, s.window = "periodic", t.window = 12)
plot(miami_tp_trends, main = "STL Decomposition of Precipitation in Miami, 1st Difference")

# stl decomposition of u10 wind
ts_miami_u10 <- ts(miami_u10, frequency = 12)
miami_u10_trends <- stl(ts_miami_u10, s.window = "periodic")
plot(miami_u10_trends, main = "Additive STL Decomposition of u10 in Miami")

# stl decomposition of v10 wind
ts_miami_v10 <- ts(miami_v10, frequency = 12)
miami_v10_trends <- stl(ts_miami_v10, s.window = "periodic")
plot(miami_v10_trends, main = "Additive STL Decomposition of v10 in Miami")
