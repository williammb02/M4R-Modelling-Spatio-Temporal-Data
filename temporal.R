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

# stl decomposition of precipitation in miami
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


# correlation analysis
# function that returns all 3 correlation measures
cors <- function(x, y){
  p <- cor(x, y, method="pearson")
  k <- cor(x, y, method="kendall")
  s <- cor(x, y, method="spearman")
  c(p,k,s)
}

# compute correlation
cors(miami_u10, miami_v10)
cors(miami_u10_rem, miami_v10_rem)

# create QF and CDF correlation plots
QFcor_plot(miami_u10, miami_v10, grid=100, xlim=c(0.05, 0.95), ylim=c(0.05, 0.95))
QFcor_plot(miami_u10_rem, miami_v10_rem, grid=100, xlim=c(0.05, 0.95), ylim=c(0.05, 0.95))
CDFcor_plot(miami_u10, miami_v10, grid=100, xlim=c(-4, 0.8), ylim=c(-2, 1.5))


# fit a distribution to u and v components
# consider u and v components before any transformations
x <- seq(0, 7, length.out=1000) 
y <- seq(-7, 5, length.out=1000)
# u component
# ghyp distribution
u_ghyp <- stepAIC.ghyp(miami_u10, silent = TRUE)
hist(u_ghyp$best.model, main="Generalised Hyperbolic, u", ylim=c(0,0.5))

# best model is essentially identical to gaussian, fit normal dist
u_norm <- fitdistr(miami_u10, "normal")
nparam <- u_norm$estimate
hist(miami_u10, freq=FALSE, breaks=50, main="Normal, u")
lines(y, dnorm(y, nparam[1], nparam[2]))

# three param weibull
u_mu <- min(miami_u10)
rescaled_u <- miami_u10 - u_mu
rescaled_u[match(0, rescaled_u)] <- 0.00000000001
rescaled_u_weib <- fitdistr(rescaled_u, "weibull")
uparam <- rescaled_u_weib$estimate
hist(miami_u10, breaks=30, freq=FALSE, main="u")
lines(y, dweibull3(y, uparam[1], uparam[2], u_mu), col="blue")
lines(y, dnorm(y, nparam[1], nparam[2]), col="red")
legend(-1, 0.4, legend=c("Normal", "3-Parameter Weibull"), 
       col = c("red", "blue"),lty=1, cex=0.8)

# v component
# try a ghyp or a three parameter weibull distribution
v_ghyp <- stepAIC.ghyp(miami_v10, silent = TRUE)
v_log <- fitdistr(miami_v10, "logistic")

# three param weibull
v_mu <- min(miami_v10)
rescaled_v <- miami_v10 - v_mu
rescaled_v[match(0, rescaled_v)] <- 0.00000000001
rescaled_v_weib <- fitdistr(rescaled_v, "weibull")
vparam <- rescaled_v_weib$estimate
hist(miami_v10, breaks=30, freq=FALSE, main="v")
lines(y, dweibull3(y, vparam[1], vparam[2], v_mu), col="blue")
lines(y, dlogis(y, v_log$estimate[1], v_log$estimate[2]), col="red")
lines(y, dghyp(y, object = v_ghyp$best.model), col="green")
legend(-4, 0.4, legend=c("Logistic", "3-Parameter Weibull", "Generalised Hyperbolic"), 
       col = c("red", "blue", "green"),lty=1, cex=0.8)

# absolute value of u and v component
# use weibull distribution
abs_u <- abs(miami_u10)
abs_v <- abs(miami_v10)

abs_u_weib <- fitdistr(abs_u, "weibull")
u_param <- abs_u_weib$estimate
hist(abs_u, breaks=30, freq=FALSE, main="2-Parameter Weibull, |u|")
lines(x, dweibull(x, u_param[1], u_param[2]), col="blue")

abs_v_weib <- fitdistr(abs_v, "weibull")
v_param <- abs_v_weib$estimate
hist(abs_v, breaks=30, freq=FALSE, main="2-Parameter Weibull, |v|")
lines(x, dweibull(x, v_param[1], v_param[2]), col="blue")

# fit a distribution to wind speed
# try ghyp
miami_speed <- sqrt(miami_u10^2 + miami_v10^2)
miami_speed_ghyp <- stepAIC.ghyp(miami_speed, silent = TRUE)
hist(miami_speed_ghyp$best.model, main="Generalised Hyperbolic, speed", ylim=c(0, 0.5))

# two parameter weibull distribution
miami_speed_weib <- fitdistr(miami_speed, "weibull")
param <- miami_speed_weib$estimate
hist(miami_speed, breaks=30, freq=FALSE, main="2-Parameter Weibull, speed")
lines(x, dweibull(x, param[1], param[2]), col="blue")

# three parameter weibull distribution
mu <- min(miami_speed)
rescaled_speed <- miami_speed - mu
rescaled_speed[match(0, rescaled_speed)] <- 0.00000000001
rescaled_weib <- fitdistr(rescaled_speed, "weibull")
param2 <- rescaled_weib$estimate

# summary plot showing distributions
hist(miami_speed, breaks=40, freq=FALSE, main="Wind Speed")
lines(x, dweibull(x, param[1], param[2]), col="red")
lines(x, dweibull3(x, param2[1], param2[2], mu), col="blue")
# lines(x, dghyp(x, object=miami_speed_ghyp$best.model), col="green")
legend(4, 0.4, legend=c("2-Parameter Weibull", "3-Parameter Weibull"),
       col = c("red", "blue"), lty=1, cex=0.8)
