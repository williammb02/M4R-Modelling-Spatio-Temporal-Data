miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)

x <- seq(0, 7, length.out=1000) 
y <- seq(-12, 10, length.out=1000)

# gaussian
u_2023_norm <- fitdistr(miami_u10_2023, "normal")
nparam <- u_2023_norm$estimate

hist(miami_u10_2023, breaks = 100, freq=FALSE)
lines(y, dnorm(y, nparam[1], nparam[2]), col="red")
lines(density(miami_u10_2023), col="blue")
legend(2, 0.14, legend=c("Normal", "Empirical"), 
       col = c("red", "blue"),lty=1, cex=0.8)
# note some more mass in the upper tail, skewed slightly compared to norm/ghyp


# gaussian mixture
v_mixture <- Mclust(miami_v10_2023)
v_mix_dens <- densityMclust(miami_v10_2023)
plot(v_mix_dens, what="density", data=miami_v10_2023, breaks = 30)
plot(v_mix_dens, what="density", data=miami_v10_2023, breaks = 100)
plot(v_mix_dens, what="BIC")
densityMclust.diagnostic(v_mix_dens, type="qq")

# parameters for model
v_weights <- v_mixture$parameters$pro
v_means <- v_mixture$parameters$mean
v_sds <- sqrt(v_mixture$parameters$variance$sigmasq)

# plot model
hist(miami_v10_2023, breaks=100, freq=FALSE, main="v")
lines(density(miami_v10_2023), col="blue")
lines(y, v_weights[1]*dnorm(y, v_means[1], v_sds[1])+v_weights[2]*dnorm(y, v_means[2], v_sds[2]), col="red")
legend(-4, 0.5, legend=c("Gaussian Mixture", "Empirical"), 
       col = c("red", "blue"),lty=1, cex=0.8)


# weibull
miami_w_weib <- fitdistr(miami_w_2023, "weibull")
param <- miami_w_weib$estimate

par(mfrow=c(1,2))
hist(miami_w_2023, breaks=100, freq=FALSE, main="Wind Speed")
lines(x, dweibull(x, param[1], param[2]), col="red")
lines(density(miami_w_2023), col="blue")
legend(3.2, 0.5, legend=c("2-Parameter Weibull", "Empirical"),
       col = c("red", "blue"), lty=1, cex=0.8)

