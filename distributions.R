miami_u10 <- import_climate_data("u10", miami_index[1], miami_index[2])
miami_v10 <- import_climate_data("v10", miami_index[1], miami_index[2])

# fit a distribution to u and v components
# consider u and v components before any transformations
x <- seq(0, 7, length.out=1000) 
y <- seq(-7, 5, length.out=1000)


# u component
# generalised hyperbolic distribution
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
hist(miami_u10, breaks=100, freq=FALSE, main="u")
lines(y, dweibull3(y, uparam[1], uparam[2], u_mu), col="blue")
lines(y, dnorm(y, nparam[1], nparam[2]), col="red")
lines(y, dghyp(y, object = u_ghyp$best.model), col="green")
lines(density(miami_u10), col="purple")
legend(-1.2, 0.4, legend=c("Normal", "3-Parameter Weibull", "Generalised Hyperbolic", "Empirical"), 
       col = c("red", "blue", "green", "purple"),lty=1, cex=0.8)

# without 3 param weib
hist(miami_u10, breaks=100, freq=FALSE, main="u")
lines(y, dnorm(y, nparam[1], nparam[2]), col="red")
lines(y, dghyp(y, object = u_ghyp$best.model), col="green")
lines(density(miami_u10), col="blue")
legend(-1.2, 0.4, legend=c("Normal", "Generalised Hyperbolic", "Empirical"), 
       col = c("red", "green", "blue"),lty=1, cex=0.8)
# note some more mass in the upper tail, skewed slightly compared to norm/ghyp

qqplot(miami_u10, qnorm(seq(0,1, by=0.005), nparam[1], nparam[2]), 
       ylab="Normal", xlab="Wind Speed", main="Normal QQ Plot")
abline(0, 1, col="purple")


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
hist(miami_v10, breaks=100, freq=FALSE, main="v")
lines(y, dweibull3(y, vparam[1], vparam[2], v_mu), col="blue")
lines(y, dlogis(y, v_log$estimate[1], v_log$estimate[2]), col="red")
lines(y, dghyp(y, object = v_ghyp$best.model), col="green")
legend(-4, 0.4, legend=c("Logistic", "3-Parameter Weibull", "Generalised Hyperbolic"), 
       col = c("red", "blue", "green"),lty=1, cex=0.8)

# try a mixture model for v
v_mixture <- Mclust(miami_v10)
v_mix_dens <- densityMclust(miami_v10)
par(mfrow=c(2,2))
plot(v_mix_dens, what="density", data=miami_v10, breaks = 30)
plot(v_mix_dens, what="density", data=miami_v10, breaks = 100)
plot(v_mix_dens, what="BIC")
densityMclust.diagnostic(v_mix_dens, type="qq")
#E: equal variance, V: variable variance

# parameters for model
v_weights <- v_mixture$parameters$pro
v_means <- v_mixture$parameters$mean
v_sds <- sqrt(v_mixture$parameters$variance$sigmasq)

# plot model
hist(miami_v10, breaks=100, freq=FALSE, main="v")
lines(density(miami_v10), col="blue")
lines(y, dlogis(y, v_log$estimate[1], v_log$estimate[2]), col="orange")
lines(y, dghyp(y, object = v_ghyp$best.model), col="purple")
lines(y, v_weights[1]*dnorm(y, v_means[1], v_sds[1])+v_weights[2]*dnorm(y, v_means[2], v_sds[2]), col="red")
legend(-4, 0.5, legend=c("Logistic", "Generalised Hyperbolic", "Gaussian Mixture", "Empirical"), 
       col = c("orange", "purple", "red", "blue"),lty=1, cex=0.8)

qqplot(miami_v10, v_weights[1]*qnorm(seq(0,1,by=0.005), v_means[1], v_sds[1])+v_weights[2]*qnorm(seq(0,1,by=0.005), v_means[2], v_sds[2]), 
       ylab="Gaussian Mixture", xlab="Wind Speed", main="Mixture Gaussian QQ Plot")
abline(0, 1, col="purple")

qqplot(miami_v10, qlogis(seq(0,1,by=0.005), v_log$estimate[1], v_log$estimate[2]), 
       ylab="Gaussian Mixture", xlab="Wind Speed", main="Mixture Gaussian QQ Plot")
abline(0, 1, col="purple")


# absolute value of u component
# use weibull distribution
abs_u <- abs(miami_u10)
abs_v <- abs(miami_v10)

abs_u_weib <- fitdistr(abs_u, "weibull")
u_param <- abs_u_weib$estimate

par(mfrow=c(1,2))
hist(abs_u, breaks=100, freq=FALSE, main="2-Parameter Weibull, |u|")
lines(x, dweibull(x, u_param[1], u_param[2]), col="red")
lines(density(abs_u), col="blue")
legend(3.2, 0.5, legend=c("Weibull", "Empirical"), 
       col=c("red", "blue"), lty=1, cex=0.8)

qqplot(abs_u, qweibull(seq(0,1,by=0.005), u_param[1], u_param[2]), 
       ylab="Weibull", xlab="Wind Speed", main="Weibull QQ Plot")
abline(0, 1, col="purple")
# might need a mixture model

# weibull mixture
abs_u_mixture <- mixfit(abs_u, family="weibull", ncomp=2)
# pi: proportion, mu: mean, sd: standard deviation
auw <- abs_u_mixture$pi
aush <- abs_u_mixture$k
ausc <- abs_u_mixture$lambda

par(mfrow=c(1,2))
hist(abs_u, breaks=100, freq=FALSE, main="Weibull Mixture, |u|")
lines(x, auw[1]*dweibull(x, aush[1], ausc[1])+auw[2]*dweibull(x, aush[2], ausc[2]), col="red")
lines(density(abs_u), col="blue")
legend(3.2, 0.5, legend=c("Weibull Mixture", "Empirical"), 
       col=c("red", "blue"), lty=1, cex=0.8)

qqplot(abs_u, auw[1]*qweibull(seq(0,1,by=0.005), aush[1], ausc[1])+auw[2]*qweibull(seq(0,1,by=0.005), aush[2], ausc[2]),
       ylab="Weibull Mixture", xlab="Wind Speed", main="Weibull Mixture QQ Plot")
abline(0, 1, col="purple")

# absolute value of v component

abs_v_weib <- fitdistr(abs_v, "weibull")
v_param <- abs_v_weib$estimate

par(mfrow=c(1,2))
hist(abs_v, breaks=100, freq=FALSE, main="2-Parameter Weibull, |v|")
lines(x, dweibull(x, v_param[1], v_param[2]), col="red")
lines(density(abs_v), col="blue")
legend(2.5, 0.8, legend=c("Weibull", "Empirical"), 
       col=c("red", "blue"), lty=1, cex=0.8)

qqplot(abs_v, qweibull(seq(0,1,by=0.005), v_param[1], v_param[2]), 
       ylab="Weibull", xlab="Wind Speed", main="Weibull QQ Plot")
abline(0, 1, col="purple")

# weibull mixture
abs_v_mixture <- mixfit(abs_v, family="weibull", ncomp=2)
# pi: proportion, mu: mean, sd: standard deviation
avw <- abs_v_mixture$pi
avsh <- abs_v_mixture$k
avsc <- abs_v_mixture$lambda

par(mfrow=c(1,2))
hist(abs_v, breaks=100, freq=FALSE, main="Weibull Mixture, |v|")
lines(x, avw[1]*dweibull(x, avsh[1], avsc[1])+avw[2]*dweibull(x, avsh[2], avsc[2]), col="red")
lines(density(abs_v), col="blue")
legend(2.5, 0.8, legend=c("Weibull", "Empirical"), 
       col=c("red", "blue"), lty=1, cex=0.8)

qqplot(abs_v, avw[1]*qweibull(seq(0,1,by=0.005), avsh[1], avsc[1])+avw[2]*qweibull(seq(0,1,by=0.005), avsh[2], avsc[2]), 
       ylab="Weibull", xlab="Wind Speed", main="Weibull QQ Plot")
abline(0, 1, col="purple")


# wind speed, norm of two components
# generalised hyperbolic distribution
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
hist(miami_speed, breaks=100, freq=FALSE, main="Wind Speed")
lines(x, dweibull(x, param[1], param[2]), col="red")
# lines(x, dweibull3(x, param2[1], param2[2], mu), col="blue")
# lines(x, dghyp(x, object=miami_speed_ghyp$best.model), col="green")
lines(density(miami_speed), col="blue")
legend(4, 0.4, legend=c("2-Parameter Weibull", "Empirical"),
       col = c("red", "blue"), lty=1, cex=0.8)

#qq plot to check the fit
qqplot(miami_speed, qweibull(seq(0,1,by=0.005), param[1], param[2]), 
       ylab="Weibull", xlab="Wind Speed", main="Weibull QQ Plot")
abline(0, 1, col="purple")


# multivariate mixture model
uv <- cbind(miami_u10, miami_v10)
uv_mix_dens <- densityMclust(uv)
plot(uv_mix_dens, what = "density", type = "hdr", data = uv, points.cex = 0.5)
uv_mixture <- Mclust(uv)
plot(uv_mix_dens, what="density", type="persp")

