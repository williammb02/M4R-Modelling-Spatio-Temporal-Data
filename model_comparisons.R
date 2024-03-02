# simulation comparison
set.seed(5)
miami_svine <- c()
for(i in 1:1000){
  miami_svine[i] <- svinecop_sim(1, 1, s_fit2)[1]
}
miami_svine <- qghyp(miami_svine, object=g11$best.model)
miami_armagarch <- ugarchsim(y_garch, n.sim=1000)
miami_armagarch <- miami_armagarch@simulation$seriesSim
miami_test <- sample(y, size=1000)

# histogram of realisations
par(mfrow=c(1,3))
hist(miami_test, breaks=30, main="Original")
hist(miami_armagarch, breaks=30, main="ARMA-GARCH")
hist(miami_svine, breaks=30, main="S-Vine")

# simulations
plot(1:1000, miami_test, type="l", col="black", xlab="", ylab="y")
lines(miami_armagarch, col="red")
lines(miami_svine, col="blue")
legend("topright", legend=c("Original", "ARMA-GARCH", "S-Vine"), col=c("black", "red", "blue"),
       lty=1, cex=0.8)

c(mean(y), mean(miami_test), mean(miami_armagarch), mean(miami_svine))
c(var(y), var(miami_test), var(miami_armagarch), var(miami_svine))



# compare moments
armagarch_means <- c()
armagarch_vars <- c()
for(i in 1:1000){
  sims <- ugarchsim(y_garch, n.sim=100)
  armagarch_means[i] <- mean(sims@simulation$seriesSim)
  armagarch_vars[i] <- var(sims@simulation$seriesSim)
}
hist(armagarch_means, breaks=100, main="Mean", xlab="Mean of Simulation")
hist(armagarch_vars, breaks=100, main="Variance", xlab="Variance of Simulation")

svine_means <- c()
svine_vars <- c()
for(i in 1:1000){
  sims <- svinecop_sim(50, 1, s_fit2)[,1]
  simsnew <- qghyp(sims, object=g11$best.model)
  svine_means[i] <- mean(simsnew)
  svine_vars[i] <- var(simsnew)
}

hist(svine_means, breaks=100, main="Mean", xlab="Mean of Simulation")
hist(svine_vars, breaks=100, main="Variance", xlab="Variance of Simulation")
