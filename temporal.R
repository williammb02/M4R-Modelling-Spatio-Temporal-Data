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
miami_tp_rem <- miami_tp_trends$time.series[,3]
par(mfrow = c(2,2))
# change between acf and pacf for an appropriate plot, pacf removes irrelevant 1.0 term
pacf(miami_tp_rem)
hist(miami_tp_rem, main = "Histogram")
plot(periodogram(as.vector(miami_tp_rem)), type = "l", ylab = "", main = "Periodogram")
qqnorm(miami_tp_rem)
qqline(miami_tp_rem, col = "red")

# fit ARIMA model to Miami total precipitation
miami_model <- auto.arima(miami_tp, ic = "aicc")
print(miami_model)
checkresiduals(miami_model)

# fit an ARIMA model to the tp remainder
miami_remainder_model <- auto.arima(miami_tp_rem, ic = "aicc")
print(miami_remainder_model)
checkresiduals(miami_remainder_model)

# fit a generalised hyperbolic distribution to the data
miami_tp_rem_ghyp <- stepAIC.ghyp(miami_tp_rem, silent = TRUE)
hist(miami_tp_rem_ghyp$best.model, main = "Histogram")

# analysis of absolute value or square of residuals
par(mfrow = c(1,2))
pacf(abs(miami_tp_rem))
pacf(miami_tp_rem^2)

# try stl decomposition of total precipitation, under first differencing
ts_miami_tp <- ts(diff(miami_tp), frequency = 12)
miami_tp_trends <- stl(ts_miami_tp, s.window = "periodic", t.window = 12)
plot(miami_tp_trends, main = "STL Decomposition of Precipitation in Miami, 1st Difference")

# stl decomposition of u10 wind
ts_miami_u10 <- ts(miami_u10, frequency = 12)
miami_u10_trends <- stl(ts_miami_u10, s.window = "periodic", t.window = length(miami_u10))
plot(miami_u10_trends, main = "Additive STL Decomposition of u10 in Miami")

# stl decomposition of v10 wind
ts_miami_v10 <- ts(miami_v10, frequency = 12)
miami_v10_trends <- stl(ts_miami_v10, s.window = "periodic", t.window = length(miami_v10))
plot(miami_v10_trends, main = "Additive STL Decomposition of v10 in Miami")

# consider all remainders
miami_tp_rem <- miami_tp_trends$time.series[, 3]
miami_u10_rem <- miami_u10_trends$time.series[, 3]
miami_v10_rem <- miami_v10_trends$time.series[, 3]

# fit ghyp distribution to remainders, and exploratory analysis
miami_urem_ghyp <- stepAIC.ghyp(miami_u10_rem, silent = TRUE)
hist(miami_urem_ghyp$best.model, ylim=c(0,0.5), main = "Histogram")
miami_vrem_ghyp <- stepAIC.ghyp(miami_v10_rem, silent = TRUE)
hist(miami_vrem_ghyp$best.model, main = "Histogram")

par(mfrow = c(2,2))
pacf(miami_u10_rem)
hist(miami_u10_rem, main = "Histogram")
plot(periodogram(as.vector(miami_u10_rem)), type = "l", ylab = "", main = "Periodogram")
qqnorm(miami_u10_rem, main = "Normal QQ Plot for u-component")
qqline(miami_u10_rem, col = "red")

par(mfrow = c(2,2))
pacf(miami_v10_rem)
hist(miami_v10_rem, main = "Histogram")
plot(periodogram(as.vector(miami_v10_rem)), type = "l", ylab = "", main = "Periodogram")
qqnorm(miami_v10_rem, main = "Normal QQ Plot for v-component")
qqline(miami_v10_rem, col = "red")

# consider features of wind decompositions
feat_stl(miami_u10, .period=12, s.window = "periodic", t.window=length(miami_u10))
feat_stl(miami_v10, .period=12, s.window = "periodic", t.window=length(miami_v10))
feat_stl(miami_tp, .period=12, s.window = "periodic", t.window=length(miami_tp))

# function that returns all 3 correlation measures
cors <- function(x, y){
  p <- cor(x, y, method="pearson")
  k <- cor(x, y, method="kendall")
  s <- cor(x, y, method="spearman")
  c(p,k,s)
}

# correlation
cors(miami_u10, miami_v10)
cors(miami_u10_rem, miami_v10_rem)

QFcor_plot(miami_u10, miami_v10, grid=100, xlim=c(0.05, 0.95), ylim=c(0.05, 0.95))
QFcor_plot(miami_u10_rem, miami_v10_rem, grid=100, xlim=c(0.05, 0.95), ylim=c(0.05, 0.95))

aa <- miami_u10
bb <- miami_v10
CDFcor_plot(aa, bb, grid=100, xlim=c(-4, 0.8), ylim=c(-2, 1.5))


