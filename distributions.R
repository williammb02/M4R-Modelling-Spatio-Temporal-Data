miami_u10 <- import_climate_data("u10", miami_index[1], miami_index[2])
miami_v10 <- import_climate_data("v10", miami_index[1], miami_index[2])

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

# empirical density f

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
lines(density(miami_u10), col="black")
legend(-1.2, 0.4, legend=c("Normal", "3-Parameter Weibull", "Generalised Hyperbolic", "Empirical"), 
       col = c("red", "blue", "green", "black"),lty=1, cex=0.8)

# without 3 param weib
hist(miami_u10, breaks=100, freq=FALSE, main="u")
lines(y, dnorm(y, nparam[1], nparam[2]), col="red")
lines(y, dghyp(y, object = u_ghyp$best.model), col="green")
lines(density(miami_u10), col="blue")
legend(-1.2, 0.4, legend=c("Normal", "Generalised Hyperbolic", "Empirical"), 
       col = c("red", "green", "blue"),lty=1, cex=0.8)


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
weights <- v_mixture$parameters$pro
means <- v_mixture$parameters$mean
sds <- v_mixture$parameters$variance$scale

mixture_density <- densityMclust(miami_v10)
hist(miami_v10, breaks=100, freq=FALSE, main="v", ylim=c(0, 0.65))
lines(density(miami_v10), col="black")
lines(y, weights[1]*dnorm(y, means[1], sds[1])+weights[2]*dnorm(y, means[2], sds[2]), col="blue")


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
