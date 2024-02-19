set.seed(5)

# generate a formula to use for finding coefficients
create_formula <- function(key_freqs) {
  sentence <- ""
  if (length(key_freqs)%%2 == 0){
    k <- length(key_freqs)/2 + 1
  } else {
    k <- round(length(key_freqs)/2, 0)
  }
  
  for (i in k:length(key_freqs)) {
    term <- paste("sin(2*pi*", key_freqs[i], "*t) + cos(2*pi*", key_freqs[i], "*t)")
    sentence <- paste(sentence, term, "+", sep = "")
  }
  return(as.formula(paste("z ~", substr(sentence, 1, nchar(sentence) - 1))))
}
# same but for monthly data
create_formulam <- function(key_freqs) {
  sentence <- ""
  if (length(key_freqs)%%2 == 0){
    k <- length(key_freqs)/2 + 1
  } else {
    k <- round(length(key_freqs)/2, 0)
  }
  
  for (i in k:length(key_freqs)) {
    term <- paste("sin(2*pi*", key_freqs[i], "*tm) + cos(2*pi*", key_freqs[i], "*tm)")
    sentence <- paste(sentence, term, "+", sep = "")
  }
  return(as.formula(paste("zm ~", substr(sentence, 1, nchar(sentence) - 1))))
}

create_formulamt <- function(key_freqs) {
  sentence <- ""
  if (length(key_freqs)%%2 == 0){
    k <- length(key_freqs)/2 + 1
  } else {
    k <- round(length(key_freqs)/2, 0)
  }
  
  for (i in k:length(key_freqs)) {
    term <- paste("sin(2*pi*", key_freqs[i], "*tm) + cos(2*pi*", key_freqs[i], "*tm)")
    sentence <- paste(sentence, term, "+", sep = "")
  }
  return(as.formula(paste("zmt ~", substr(sentence, 1, nchar(sentence) - 1))))
}

t <- 1:length(miami_w_2023)
freq <- seq(from = -0.5, to = 0.5, length = length(miami_w_2023))
tm <- 1:length(miami_w)
freqm <- seq(from = -0.5, to = 0.5, length = length(miami_w))

# fit a personal stl model to miami wind speed
# aim for an additive model

# fit the trend
# acf(miami_w_2023)

# hourly wind speed decompositions
quadfit <- lm(miami_w_2023 ~ poly(t, 2, raw=TRUE))
z <- miami_w_2023 - quadfit$fitted.values
key_freqs <- c()
for(i in 1:length(my_periodogram(z))){
  if(my_periodogram(z)[i] > 170){
    key_freqs <- c(key_freqs, freq[i])
  }
}
key_freqs <- key_freqs[11:20]
seasonqfit <- lm(create_formula(key_freqs))
y <- z - seasonqfit$fitted.values
# fit arima model to y
y_model <- auto.arima(y, ic = "aicc")
# fit arma garch model to y 
yspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                    mean.model=list(armaOrder = c(1, 4), include.mean=TRUE))
y_garch <- ugarchfit(spec = yspec, data = y, solver="lbfgs")
y_garch@fit$solver$sol$par
y_gar_res <- y_garch@fit$residuals
y_armac <- y_garch@fit$coef[2:6]
y_garchc <- y_garch@fit$coef[7:9]
# tampa wind speed
# quadratic
quadfit2 <- lm(tampa_w_2023 ~ poly(t, 2, raw=TRUE))
z2 <- tampa_w_2023 - quadfit2$fitted.values
key_freqs2 <- c()
for(i in 1:length(my_periodogram(z2))){
  if(my_periodogram(z2)[i] > 125){
    key_freqs2 <- c(key_freqs2, freq[i])
  }
}
key_freqs2 <- key_freqs2[11:20]
seasonqfit2 <- lm(create_formula(key_freqs2))
y2 <- z2 - seasonqfit2$fitted.values
# fit arima model to y2
y2_model <- auto.arima(y2, ic = "aicc")

y2spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(5, 2), include.mean=TRUE))
y2_garch <- ugarchfit(spec = y2spec, data = y2, solver="hybrid")
y2_garch@fit$solver$sol$par
y2_armac <- y2_garch@fit$solver$sol$par[2:8]
y2_garchc <- y2_garch@fit$solver$sol$par[9:11]
y2_gar_res <- y2_garch@fit$residuals
# tallahassee ARMA-GARCH model
quadfit3 <- lm(talla_w_2023 ~ poly(t, 2, raw=TRUE))
z3 <- talla_w_2023 - quadfit3$fitted.values
key_freqs3 <- c()
for(i in 1:length(my_periodogram(z3))){
  if(my_periodogram(z3)[i] > 31.5){
    key_freqs3 <- c(key_freqs3, freq[i])
  }
}
key_freqs3 <- key_freqs3[11:20]
seasonqfit3 <- lm(create_formula(key_freqs3))
y3 <- z3 - seasonqfit3$fitted.values
# fit arima model to y3
y3_model <- auto.arima(y3, ic = "aicc")
# fit arma garch model to y3 using arima choice 
y3spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(1,1), include.mean=TRUE))
y3_garch <- ugarchfit(spec = y3spec, data = y3, solver="hybrid")
y3_garch@fit$solver$sol$par
y3_armac <- y3_garch@fit$solver$sol$par[2:3]
y3_garchc <- y3_garch@fit$solver$sol$par[4:6]
y3_gar_res <- y3_garch@fit$residuals
# jackson
quadfit4 <- lm(jack_w_2023 ~ poly(t, 2, raw=TRUE))
z4 <- jack_w_2023 - quadfit4$fitted.values
key_freqs4 <- c()
for(i in 1:length(my_periodogram(z4))){
  if(my_periodogram(z4)[i] > 47){
    key_freqs4 <- c(key_freqs4, freq[i])
  }
}
key_freqs4 <- key_freqs4[11:20]
seasonqfit4 <- lm(create_formula(key_freqs4))
y4 <- z4 - seasonqfit4$fitted.values
# fit arima model to y4
y4_model <- auto.arima(y4, ic = "aicc")
# fit arma garch model to y4 using arima choice 
y4spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(2,1), include.mean=TRUE))
y4_garch <- ugarchfit(spec = y4spec, data = y4, solver="hybrid")
y4_garch@fit$solver$sol$par
y4_gar_res <- y4_garch@fit$residuals

# orlando
quadfit5 <- lm(orla_w_2023 ~ poly(t, 2, raw=TRUE))
z5 <- orla_w_2023 - quadfit5$fitted.values
key_freqs5 <- c()
for(i in 1:length(my_periodogram(z5))){
  if(my_periodogram(z5)[i] > 77){
    key_freqs5 <- c(key_freqs5, freq[i])
  }
}
key_freqs5 <- key_freqs5[11:20]
seasonqfit5 <- lm(create_formula(key_freqs5))
y5 <- z5 - seasonqfit5$fitted.values
# fit arima model to y5
y5_model <- auto.arima(y5, ic = "aicc")
# fit arma garch model to y5 using arima choice 
y5spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(1,1), include.mean=TRUE))
y5_garch <- ugarchfit(spec = y5spec, data = y5, solver="hybrid")
y5_garch@fit$solver$sol$par
y5_gar_res <- y5_garch@fit$residuals

# fort meyers
quadfit6 <- lm(fort_w_2023 ~ poly(t, 2, raw=TRUE))
z6 <- fort_w_2023 - quadfit6$fitted.values
key_freqs6 <- c()
for(i in 1:length(my_periodogram(z6))){
  if(my_periodogram(z6)[i] > 97.75){
    key_freqs6 <- c(key_freqs6, freq[i])
  }
}
key_freqs6 <- key_freqs6[11:20]
seasonqfit6 <- lm(create_formula(key_freqs6))
y6 <- z6 - seasonqfit6$fitted.values
# fit arima model to y6
y6_model <- auto.arima(y6, ic = "aicc")
# fit arma garch model to y6 using arima choice 
y6spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(3,1), include.mean=TRUE))
y6_garch <- ugarchfit(spec = y6spec, data = y6, solver="hybrid")
y6_garch@fit$solver$sol$par
y6_gar_res <- y6_garch@fit$residuals

# hourly precipitation
# miami
quadfitt <- lm(miami_tp_2023 ~ poly(t, 2, raw=TRUE))
zt <- miami_tp_2023 - quadfitt$fitted.values
key_freqst <- c()
for(i in 1:length(my_periodogram(zt))){
  if(my_periodogram(zt)[i] > 0.0004){
    key_freqst <- c(key_freqst, freq[i])
  }
}
key_freqst <- key_freqst[11:20]
seasonqfitt <- lm(create_formula(key_freqst))
yt <- zt - seasonqfitt$fitted.values

# tampa
quadfit2t <- lm(tampa_tp_2023 ~ poly(t, 2, raw=TRUE))
zt2 <- tampa_tp_2023 - quadfit2t$fitted.values
key_freqs2t <- c()
for(i in 1:length(my_periodogram(zt2))){
  if(my_periodogram(zt2)[i] > 0.0002){
    key_freqs2t <- c(key_freqs2t, freq[i])
  }
}
key_freqs2t <- key_freqs2t[13:22]
seasonqfit2t <- lm(create_formula(key_freqs2t))
yt2 <- zt2 - seasonqfit2t$fitted.values

# tallahassee
quadfit3t <- lm(talla_tp_2023 ~ poly(t, 2, raw=TRUE))
zt3 <- talla_tp_2023 - quadfit3t$fitted.values
key_freqs3t <- c()
for(i in 1:length(my_periodogram(zt3))){
  if(my_periodogram(zt2)[i] > 0.0002){
    key_freqs3t <- c(key_freqs3t, freq[i])
  }
}
key_freqs3t <- key_freqs3t[13:22]
seasonqfit3t <- lm(create_formula(key_freqs3t))
yt3 <- zt3 - seasonqfit3t$fitted.values

#jacksonville
quadfit4t <- lm(jack_tp_2023 ~ poly(t, 2, raw=TRUE))
zt4 <- jack_tp_2023 - quadfit4t$fitted.values
key_freqs4t <- c()
for(i in 1:length(my_periodogram(zt4))){
  if(my_periodogram(zt4)[i] > 0.0003){
    key_freqs4t <- c(key_freqs4t, freq[i])
  }
}
key_freqs4t <- key_freqs4t[11:20]
seasonqfit4t <- lm(create_formula(key_freqs4t))
yt4 <- zt4 - seasonqfit4t$fitted.values

# orlando
quadfit5t <- lm(orla_tp_2023 ~ poly(t, 2, raw=TRUE))
zt5 <- orla_tp_2023 - quadfit5t$fitted.values
key_freqs5t <- c()
for(i in 1:length(my_periodogram(zt5))){
  if(my_periodogram(zt5)[i] > 0.00015){
    key_freqs5t <- c(key_freqs5t, freq[i])
  }
}
key_freqs5t <- key_freqs5t[19:28]
seasonqfit5t <- lm(create_formula(key_freqs5t))
yt5 <- zt5 - seasonqfit5t$fitted.values

# fort myers
quadfit6t <- lm(fort_tp_2023 ~ poly(t, 2, raw=TRUE))
zt6 <- fort_tp_2023 - quadfit6t$fitted.values
key_freqs6t <- c()
for(i in 1:length(my_periodogram(zt6))){
  if(my_periodogram(zt6)[i] > 0.0001){
    key_freqs6t <- c(key_freqs6t, freq[i])
  }
}
key_freqs6t <- key_freqs6t[21:30]
seasonqfit6t <- lm(create_formula(key_freqs6t))
yt6 <- zt6 - seasonqfit6t$fitted.values




# monthly approach
# miami
quadfitm <- lm(miami_w ~ poly(tm, 2, raw=TRUE))
zm <- miami_w - quadfitm$fitted.values
key_freqsm <- c()
for(i in 1:length(my_periodogram(zm))){
  if(my_periodogram(zm)[i] > 3){
    key_freqsm <- c(key_freqsm, freqm[i])
  }
}
key_freqsm <- key_freqsm[11:20]
seasonqfitm <- lm(create_formulam(key_freqsm))
ym <- zm - seasonqfitm$fitted.values
# fit arima model to y
y_modelm <- auto.arima(ym, ic = "aicc")
# fit arma garch model to y 
yspecm <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                    mean.model=list(armaOrder = c(2, 1), include.mean=TRUE))
y_garchm <- ugarchfit(spec = yspecm, data = ym, solver="hybrid")
y_gar_resm <- y_garchm@fit$residuals
# tampa
quadfitm2 <- lm(tampa_w ~ poly(tm, 2, raw=TRUE))
zm2 <- tampa_w - quadfitm2$fitted.values
key_freqsm2 <- c()
for(i in 1:length(my_periodogram(zm2))){
  if(my_periodogram(zm2)[i] > 1.25){
    key_freqsm2 <- c(key_freqsm2, freqm[i])
  }
}
key_freqsm2 <- key_freqsm2[11:20]
seasonqfitm2 <- lm(create_formulam(key_freqsm2))
ym2 <- zm2 - seasonqfitm2$fitted.values
# fit arima model to y
y_modelm2 <- auto.arima(ym2, ic = "aicc")
# fit arma garch model to y 
yspecm2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                     mean.model=list(armaOrder = c(2, 1), include.mean=TRUE))
y_garchm2 <- ugarchfit(spec = yspecm2, data = ym2, solver="hybrid")
y_gar_resm2 <- y_garchm2@fit$residuals
# tallahassee
quadfitm3 <- lm(tallahassee_w ~ poly(tm, 2, raw=TRUE))
zm3 <- tallahassee_w - quadfitm3$fitted.values
key_freqsm3 <- c()
for(i in 1:length(my_periodogram(zm3))){
  if(my_periodogram(zm3)[i] > 0.25){
    key_freqsm3 <- c(key_freqsm3, freqm[i])
  }
}
key_freqsm3 <- key_freqsm3[23:32]
seasonqfitm3 <- lm(create_formulam(key_freqsm3))
ym3 <- zm3 - seasonqfitm3$fitted.values
# fit arima model to y
y_modelm3 <- auto.arima(ym3, ic = "aicc")
# fit arma garch model to y 
yspecm3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model=list(armaOrder = c(1, 4), include.mean=TRUE))
y_garchm3 <- ugarchfit(spec = yspecm3, data = ym3, solver="hybrid")
y_gar_resm3 <- y_garchm3@fit$residuals
# jacksonville
quadfitm4 <- lm(jackson_w ~ poly(tm, 2, raw=TRUE))
zm4 <- jackson_w - quadfitm4$fitted.values
key_freqsm4 <- c()
for(i in 1:length(my_periodogram(zm4))){
  if(my_periodogram(zm4)[i] > 0.5){
    key_freqsm4 <- c(key_freqsm4, freqm[i])
  }
}
key_freqsm4 <- key_freqsm4[17:26]
seasonqfitm4 <- lm(create_formulam(key_freqsm4))
ym4 <- zm4 - seasonqfitm4$fitted.values
# fit arima model to y
y_modelm4 <- auto.arima(ym4, ic = "aicc")
# fit arma garch model to y 
yspecm4 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model=list(armaOrder = c(5, 1), include.mean=TRUE))
y_garchm4 <- ugarchfit(spec = yspecm4, data = ym4, solver="hybrid")
y_gar_resm4 <- y_garchm4@fit$residuals
# orlando
quadfitm5 <- lm(orlando_w ~ poly(tm, 2, raw=TRUE))
zm5 <- orlando_w - quadfitm5$fitted.values
key_freqsm5 <- c()
for(i in 1:length(my_periodogram(zm5))){
  if(my_periodogram(zm5)[i] > 0.9){
    key_freqsm5 <- c(key_freqsm5, freqm[i])
  }
}
key_freqsm5 <- key_freqsm5[13:22]
seasonqfitm5 <- lm(create_formulam(key_freqsm5))
ym5 <- zm5 - seasonqfitm5$fitted.values
# fit arima model to y
y_modelm5 <- auto.arima(ym5, ic = "aicc")
# fit arma garch model to y 
yspecm5 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model=list(armaOrder = c(2, 3), include.mean=TRUE))
y_garchm5 <- ugarchfit(spec = yspecm5, data = ym5, solver="hybrid")
y_gar_resm5 <- y_garchm5@fit$residuals
# fort myers
quadfitm6 <- lm(fort_w ~ poly(tm, 2, raw=TRUE))
zm6 <- fort_w - quadfitm6$fitted.values
key_freqsm6 <- c()
for(i in 1:length(my_periodogram(zm6))){
  if(my_periodogram(zm6)[i] > 1){
    key_freqsm6 <- c(key_freqsm6, freqm[i])
  }
}
key_freqsm6 <- key_freqsm6[19:28]
seasonqfitm6 <- lm(create_formulam(key_freqsm6))
ym6 <- zm6 - seasonqfitm6$fitted.values
# fit arima model to y
y_modelm6 <- auto.arima(ym6, ic = "aicc")
# fit arma garch model to y 
yspecm6 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model=list(armaOrder = c(3, 1), include.mean=TRUE))
y_garchm6 <- ugarchfit(spec = yspecm6, data = ym6, solver="hybrid")
y_gar_resm6 <- y_garchm6@fit$residuals

# monthly precipitation
# miami
quadfitmt <- lm(miami_tp ~ poly(tm, 2, raw=TRUE))
zmt <- miami_tp - quadfitmt$fitted.values
key_freqsmt <- c()
for(i in 1:length(my_periodogram(zmt))){
  if(my_periodogram(zmt)[i] > 0.00001){
    key_freqsmt <- c(key_freqsmt, freqm[i])
  }
}
key_freqsmt <- key_freqsmt[25:34]
seasonqfitmt <- lm(create_formulamt(key_freqsmt))
ymt <- zmt - seasonqfitmt$fitted.values

# tampa
quadfitmt2 <- lm(tampa_tp ~ poly(tm, 2, raw=TRUE))
zmt2 <- tampa_tp - quadfitmt2$fitted.values
key_freqsmt2 <- c()
for(i in 1:length(my_periodogram(zmt2))){
  if(my_periodogram(zmt2)[i] > 0.00001){
    key_freqsmt2 <- c(key_freqsmt2, freqm[i])
  }
}
key_freqsmt2 <- key_freqsmt2[23:32]
seasonqfitmt2 <- lm(create_formulamt(key_freqsmt2))
ymt2 <- zmt2 - seasonqfitmt2$fitted.values

# tallahassee
quadfitmt3 <- lm(tallahassee_tp ~ poly(tm, 2, raw=TRUE))
zmt3 <- tallahassee_tp - quadfitmt3$fitted.values
key_freqsmt3 <- c()
for(i in 1:length(my_periodogram(zmt3))){
  if(my_periodogram(zmt3)[i] > 0.00001){
    key_freqsmt3 <- c(key_freqsmt3, freqm[i])
  }
}
key_freqsmt3 <- key_freqsmt3[17:26]
seasonqfitmt3 <- lm(create_formulamt(key_freqsmt3))
ymt3 <- zmt3 - seasonqfitmt3$fitted.values

# jacksonville
quadfitmt4 <- lm(jackson_tp ~ poly(tm, 2, raw=TRUE))
zmt4 <- jackson_tp - quadfitmt4$fitted.values
key_freqsmt4 <- c()
for(i in 1:length(my_periodogram(zmt4))){
  if(my_periodogram(zmt4)[i] > 0.00001){
    key_freqsmt4 <- c(key_freqsmt4, freqm[i])
  }
}
key_freqsmt4 <- key_freqsmt4[12:21]
seasonqfitmt4 <- lm(create_formulamt(key_freqsmt4))
ymt4 <- zmt4 - seasonqfitmt4$fitted.values

# orlando
quadfitmt5 <- lm(orlando_tp ~ poly(tm, 2, raw=TRUE))
zmt5 <- orlando_tp - quadfitmt5$fitted.values
key_freqsmt5 <- c()
for(i in 1:length(my_periodogram(zmt5))){
  if(my_periodogram(zmt5)[i] > 0.00001){
    key_freqsmt5 <- c(key_freqsmt5, freqm[i])
  }
}
key_freqsmt5 <- key_freqsmt5[18:27]
seasonqfitmt5 <- lm(create_formulamt(key_freqsmt5))
ymt5 <- zmt5 - seasonqfitmt5$fitted.values

# fort myers
quadfitmt6 <- lm(fort_tp ~ poly(tm, 2, raw=TRUE))
zmt6 <- fort_tp - quadfitmt6$fitted.values
key_freqsmt6 <- c()
for(i in 1:length(my_periodogram(zmt6))){
  if(my_periodogram(zmt6)[i] > 0.00001){
    key_freqsmt6 <- c(key_freqsmt6, freqm[i])
  }
}
key_freqsmt6 <- key_freqsmt6[21:30]
seasonqfitmt6 <- lm(create_formulamt(key_freqsmt6))
ymt6 <- zmt6 - seasonqfitmt6$fitted.values



# plots
plot(miami_w_2023, type="l", main="Quadratic Trend")
lines(quadfit$fitted.values, col="red")

plot(freq, my_periodogram(z), type="l")

plot(z, type="l", main="Quadratic Trend with Seasonal Component")
lines(seasonqfit$fitted.values, col="red")

par(mfrow = c(2,2))
plot(y, type="l", main="Remainder with Quadratic Trend")
hist(y, breaks=100, freq=FALSE, main="Histogram")
lines(density(y), col="red")
qqnorm(y, main="QQ Plot")
acf(y)

# plots
plot(tampa_w_2023, type="l", main="Quadratic Trend")
lines(quadfit2$fitted.values, col="red")

plot(freq, my_periodogram(z2), type="l")

plot(z2, type="l", main="Quadratic Trend with Seasonal Component")
lines(seasonqfit2$fitted.values, col="red")

par(mfrow = c(2,2))
plot(y2, type="l", main="Remainder with Quadratic Trend")
hist(y2, breaks=100, freq=FALSE, main="Histogram")
lines(density(y2), col="red")
qqnorm(y2, main="QQ Plot")
acf(y2)

# ARIMA RESIDUALS
# fit mvn distribution to residuals
res_mvn <- mvn("XXX", cbind(y_model$residuals, y2_model$residuals))

# fit ghyp distribution to residuals
res_ghyp <- stepAIC.ghyp(cbind(y_model$residuals, y2_model$residuals), silent = TRUE)

res1_ghyp <- stepAIC.ghyp(y_model$residuals, silent=TRUE)
res2_ghyp <- stepAIC.ghyp(y2_model$residuals, silent=TRUE)

# ljung box test
checkresiduals(y_model)
checkresiduals(y2_model)


# ARIMA forecasting
t <- 5832
y11 <- y[t] + y_model$coef[1]*y_model$residuals[t] + y_model$coef[2]*y_model$residuals[t-1] + y_model$coef[3]*y_model$residuals[t-2] + y_model$coef[4]*y_model$residuals[t-3]
y12 <- y11 + y_model$coef[2]*y_model$residuals[t] + y_model$coef[3]*y_model$residuals[t-1] + y_model$coef[4]*y_model$residuals[t-2]
y13 <- y12 + y_model$coef[3]*y_model$residuals[t] + y_model$coef[4]*y_model$residuals[t-1]
y14 <- y13 + y_model$coef[4]*y_model$residuals[t] 
y15 <- y14
y21 <- y2_model$coef[1]*y2[t] + y2_model$coef[2]*y2[t-1] + y2_model$coef[3]*y2[t-2] + y2_model$coef[4]*y2[t-3] + y2_model$coef[5]*y2[t-4] + y2_model$coef[6]*y2_model$residuals[t] + y2_model$coef[7]*y2_model$residuals[t-1]
y22 <- y2_model$coef[1]*y21 + y2_model$coef[2]*y2[t] + y2_model$coef[3]*y2[t-1] + y2_model$coef[4]*y2[t-2] + y2_model$coef[5]*y2[t-3] + y2_model$coef[7]*y2_model$residuals[t]
y23 <- y2_model$coef[1]*y22 + y2_model$coef[2]*y21 + y2_model$coef[3]*y2[t] + y2_model$coef[4]*y2[t-1] + y2_model$coef[5]*y2[t-2]
y24 <- y2_model$coef[1]*y23 + y2_model$coef[2]*y22 + y2_model$coef[3]*y21 + y2_model$coef[4]*y2[t] + y2_model$coef[5]*y2[t-1]
y25 <- y2_model$coef[1]*y24 + y2_model$coef[2]*y23 + y2_model$coef[3]*y22 + y2_model$coef[4]*y21 + y2_model$coef[5]*y2[t]

# ARIMA check forecasting inside time series lines up with what was observed
# consider after the first 1000 observations
k <- 1000
y11 <- y[k] + y_model$coef[1]*y_model$residuals[k] + y_model$coef[2]*y_model$residuals[k-1] + y_model$coef[3]*y_model$residuals[k-2] + y_model$coef[4]*y_model$residuals[k-3]
y12 <- y11 + y_model$coef[2]*y_model$residuals[k] + y_model$coef[3]*y_model$residuals[k-1] + y_model$coef[4]*y_model$residuals[k-2]
y13 <- y12 + y_model$coef[3]*y_model$residuals[k] + y_model$coef[4]*y_model$residuals[k-1]
y14 <- y13 + y_model$coef[4]*y_model$residuals[k] 
y15 <- y14
y21 <- y2_model$coef[1]*y2[k] + y2_model$coef[2]*y2[k-1] + y2_model$coef[3]*y2[k-2] + y2_model$coef[4]*y2[k-3] + y2_model$coef[5]*y2[k-4] + y2_model$coef[6]*y2_model$residuals[k] + y2_model$coef[7]*y2_model$residuals[k-1]
y22 <- y2_model$coef[1]*y21 + y2_model$coef[2]*y2[k] + y2_model$coef[3]*y2[k-1] + y2_model$coef[4]*y2[k-2] + y2_model$coef[5]*y2[k-3] + y2_model$coef[7]*y2_model$residuals[k]
y23 <- y2_model$coef[1]*y22 + y2_model$coef[2]*y21 + y2_model$coef[3]*y2[k] + y2_model$coef[4]*y2[k-1] + y2_model$coef[5]*y2[k-2]
y24 <- y2_model$coef[1]*y23 + y2_model$coef[2]*y22 + y2_model$coef[3]*y21 + y2_model$coef[4]*y2[k] + y2_model$coef[5]*y2[k-1]
y25 <- y2_model$coef[1]*y24 + y2_model$coef[2]*y23 + y2_model$coef[3]*y22 + y2_model$coef[4]*y21 + y2_model$coef[5]*y2[k]

y1preds <- rep(y15, times=495)
y1preds <- c(y11, y12, y13, y14, y15, y1preds)

y2preds <- c(y21, y22, y23, y24, y25, rep(0, times=495))
for(i in 6:500){
  y2preds[i] <- y2_model$coef[1]*y2preds[i-1] + y2_model$coef[2]*y2preds[i-2] + y2_model$coef[3]*y2preds[i-3] + y2_model$coef[4]*y2preds[i-4] + y2_model$coef[5]*y2preds[i-5]
}

yts <- 1:1500
plot(yts, y[1:1500], type="l", xlab="Index", ylab="Speed", main="Miami")
lines(yts[1001:1500], y1preds, col="blue")
abline(v = 1000, col="red", lty=2)

plot(yts, y2[1:1500], type="l", xlab="Index", ylab="Speed", main="Tampa")
lines(yts[1001:1500], y2preds, col="blue")
abline(v = 1000, col="red", lty=2)


# ARMA GARCH

# ARMA GARCH residuals distribution
garch_res_ghyp <- stepAIC.ghyp(cbind(y_garch@fit$residuals, y2_garch@fit$residuals), silent = TRUE)
garch_res_ghyp$best.model

# ARMA GARCH forecasting
t <- 5832

y11g <- y_armac[1]*y[t] + y_armac[2]*y_gar_res[t] + y_armac[3]*y_gar_res[t-1] + y_armac[4]*y_gar_res[t-2] + y_armac[5]*y_gar_res[t-3]
y12g <- y_armac[1]*y11 + y_armac[3]*y_gar_res[t] + y_armac[4]*y_gar_res[t-1] + y_armac[5]*y_gar_res[t-2]
y13g <- y_armac[1]*y12 + y_armac[4]*y_gar_res[t] + y_armac[5]*y_gar_res[t-1]
y14g <- y_armac[1]*y13 + y_armac[5]*y_gar_res[t] 
y15g <- y_armac[1]*y14
y21g <- y2_armac[1]*y2[t] + y2_armac[2]*y2[t-1] + y2_armac[3]*y2[t-2] + y2_armac[4]*y2[t-3] + y2_armac[5]*y2[t-4] + y2_armac[6]*y2_gar_res[t] + y2_armac[7]*y2_gar_res[t-1]
y22g <- y2_armac[1]*y21 + y2_armac[2]*y2[t] + y2_armac[3]*y2[t-1] + y2_armac[4]*y2[t-2] + y2_armac[5]*y2[t-3] + y2_armac[7]*y2_gar_res[t]
y23g <- y2_armac[1]*y22 + y2_armac[2]*y21 + y2_armac[3]*y2[t] + y2_armac[4]*y2[t-1] + y2_armac[5]*y2[t-2]
y24g <- y2_armac[1]*y23 + y2_armac[2]*y22 + y2_armac[3]*y21 + y2_armac[4]*y2[t] + y2_armac[5]*y2[t-1]
y25g <- y2_armac[1]*y24 + y2_armac[2]*y23 + y2_armac[3]*y22 + y2_armac[4]*y21 + y2_armac[5]*y2[t]

y1predsg <- c(y11g, y12g, y13g, y14g, y15g, rep(0, times=495))
for(i in 6:500){
  y1predsg[i] <- y_armac[1]*y1predsg[i-1]
}

y2predsg <- c(y21g, y22g, y23g, y24g, y25g, rep(0, times=495))
for(i in 6:500){
  y2predsg[i] <- y2_armac[1]*y2preds[i-1] + y2_armac[2]*y2preds[i-2] + y2_armac[3]*y2preds[i-3] + y2_armac[4]*y2preds[i-4] + y2_armac[5]*y2preds[i-5]
}

yts <- 1:1500
plot(yts, y[1:1500], type="l", xlab="Index", ylab="Speed", main="Miami")
lines(yts[1001:1500], y1predsg, col="blue")
abline(v = 1000, col="red", lty=2)

plot(yts, y2[1:1500], type="l", xlab="Index", ylab="Speed", main="Tampa")
lines(yts[1001:1500], y2preds, col="blue")
abline(v = 1000, col="red", lty=2)



#ARMA GARCH conditional variance forecasting
sigma1 <- coredata(sigma(y_garch))[t,1]
sigma2 <- coredata(sigma(y2_garch))[t,1]

s11 <- y_garchc[1] + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
s12 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3])) + (y_garchc[2] + y_garchc[3])*(y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
s13 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^2
s14 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2 + (y_garchc[2] + y_garchc[3])^3) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^3
s15 <- y_garchc[1]*(1 + (y_garchc[2] + y_garchc[3]) + (y_garchc[2] + y_garchc[3])^2 + (y_garchc[2] + y_garchc[3])^3 + (y_garchc[2] + y_garchc[3])^4) + (y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)*(y_garchc[2] + y_garchc[3])^4
s21 <- y2_garchc[1] + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
s22 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3])) + (y2_garchc[2] + y2_garchc[3])*(y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
s23 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^2
s24 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2 + (y2_garchc[2] + y2_garchc[3])^3) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^3
s25 <- y2_garchc[1]*(1 + (y2_garchc[2] + y2_garchc[3]) + (y2_garchc[2] + y2_garchc[3])^2 + (y2_garchc[2] + y2_garchc[3])^3 + (y2_garchc[2] + y2_garchc[3])^4) + (y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)*(y2_garchc[2] + y2_garchc[3])^4

# check forecasting inside time series lines up with what was observed
k <- 1000

y11g <- y_armac[1]*y[k] + y_armac[2]*y_gar_res[k] + y_armac[3]*y_gar_res[k-1] + y_armac[4]*y_gar_res[k-2] + y_armac[5]*y_gar_res[k-3]
y12g <- y_armac[1]*y11 + y_armac[3]*y_gar_res[k] + y_armac[4]*y_gar_res[k-1] + y_armac[5]*y_gar_res[k-2]
y13g <- y_armac[1]*y12 + y_armac[4]*y_gar_res[k] + y_armac[5]*y_gar_res[k-1]
y14g <- y_armac[1]*y13 + y_armac[5]*y_gar_res[k] 
y15g <- y_armac[1]*y14
y21g <- y2_armac[1]*y2[k] + y2_armac[2]*y2[k-1] + y2_armac[3]*y2[k-2] + y2_armac[4]*y2[k-3] + y2_armac[5]*y2[k-4] + y2_armac[6]*y2_gar_res[k] + y2_armac[7]*y2_gar_res[k-1]
y22g <- y2_armac[1]*y21 + y2_armac[2]*y2[k] + y2_armac[3]*y2[k-1] + y2_armac[4]*y2[k-2] + y2_armac[5]*y2[k-3] + y2_armac[7]*y2_gar_res[k]
y23g <- y2_armac[1]*y22 + y2_armac[2]*y21 + y2_armac[3]*y2[k] + y2_armac[4]*y2[k-1] + y2_armac[5]*y2[k-2]
y24g <- y2_armac[1]*y23 + y2_armac[2]*y22 + y2_armac[3]*y21 + y2_armac[4]*y2[k] + y2_armac[5]*y2[k-1]
y25g <- y2_armac[1]*y24 + y2_armac[2]*y23 + y2_armac[3]*y22 + y2_armac[4]*y21 + y2_armac[5]*y2[k]


# conditional variance formula
miami_s_fc <- function(tau){
  t <- 5832
  summand <- rep(0, times = tau)
  summand[1] <- 1
  if(tau > 1){
    for(i in 2:tau){
      summand[i] <- (y_garchc[2] + y_garchc[3])^(i-1)
    }
  }
  sig <- y_garchc[1]*sum(summand) + (y_garchc[2] + y_garchc[3])^(tau-1)*(y_garchc[2]*y_gar_res[t]^2 + y_garchc[3]*sigma1^2)
  return(sig)
}

tampa_s_fc <- function(tau){
  t <- 5832
  summand <- rep(0, times = tau)
  summand[1] <- 1
  if(tau > 1){
    for(i in 2:tau){
      summand[i] <- (y2_garchc[2] + y2_garchc[3])^(i-1)
    }
  }
  sig <- y2_garchc[1]*sum(summand) + (y2_garchc[2] + y2_garchc[3])^(tau-1)*(y2_garchc[2]*y2_gar_res[t]^2 + y2_garchc[3]*sigma2^2)
  return(sig)
}

xtss <- 1:(length(miami_w_2023)+48)
miami_s <- rep(0, times = 48)
tampa_s <- rep(0, times = 48)
for(i in 1:48){
  miami_s[i] <- miami_s_fc(i)
  tampa_s[i] <- tampa_s_fc(i)
}

# plot showing last 1 week and then prediction for next week
plot(xtss[(t-168):length(xtss)], c(coredata(sigma(y_garch))[(t-168):t,1], miami_s), 
     type="l", xlab="Index", ylab="Variance", main="Miami")
abline(v = length(miami_w_2023), col="red", lty=2)

plot(xtss[(t-168):length(xtss)], c(coredata(sigma(y2_garch))[(t-168):t,1], tampa_s), 
     type="l", xlab="Index", ylab="Variance", main="Tampa")
abline(v = length(tampa_w_2023), col="red", lty=2)





# VAR model
ys_var <- VAR(cbind(y, y2), p=1, type="none")
# VAR prediction
phi_var <- matrix(c(ys_var$varresult$y$coefficients, ys_var$varresult$y2$coefficients),
                     nrow=2, ncol=2, byrow=TRUE)
Y_fin <- c(y[5832], y2[5832])

var_forecast <- function(tau){
  fore <- Y_fin
  for (i in 1:tau){
    fore <- phi_var %*% fore
  }
  return(fore)
}

# create a plot showing the time series and then the forecasts
miamiforecast <- c(rep(0, times = 504))
tampaforecast <- c(rep(0, times = 504))

xts <- 1:(length(miami_w_2023)+504)

for(i in 1:504){
  miamiforecast[i] <- var_forecast(i)[1]
  tampaforecast[i] <- var_forecast(i)[2]
}

par(mfrow=c(1,2))
plot(xts[5000:length(xts)], c(y, miamiforecast)[5000:length(xts)], 
     type="l", main="Miami", xlab="Index", ylab="Speed")
abline(v = length(miami_w_2023), col="red", lty=2)

plot(xts[5000:length(xts)], c(y2, tampaforecast)[5000:length(xts)], 
     type="l", main="Tampa", xlab="Index", ylab="Speed")
abline(v = length(tampa_w_2023), col="red", lty=2)

y_var_res <- ys_var$varresult$y$residuals
y2_var_res <- ys_var$varresult$y2$residuals

var_res_mvn <- mvn("XXX", cbind(y_var_res, y2_var_res))
var_res_ghyp <- stepAIC.ghyp(cbind(y_var_res, y2_var_res), silent=TRUE)


# miami, talla
garch2_res_ghyp <- stepAIC.ghyp(cbind(y_garch@fit$residuals, y3_garch@fit$residuals), silent = TRUE)
garch2_res_ghyp$best.model

# tampa, talla
garch3_res_ghyp <- stepAIC.ghyp(cbind(y2_garch@fit$residuals, y3_garch@fit$residuals), silent = TRUE)
garch3_res_ghyp$best.model


# try fitting another seasonal component to the residuals
key_freqs11 <- c()
for(i in 1:length(my_periodogram(y_gar_res))){
  if(my_periodogram(y_gar_res)[i] > 1.9){
    key_freqs11 <- c(key_freqs11, freq[i])
  }
}
key_freqs11 <- key_freqs11[12:21]
seasongar1fit <- lm(create_formula(key_freqs11))
y_gar_res_new <- y_gar_res - seasongar1fit$fitted.values

key_freqs22 <- c()
for(i in 1:length(my_periodogram(y2_gar_res))){
  if(my_periodogram(y2_gar_res)[i] > 1.8){
    key_freqs22 <- c(key_freqs22, freq[i])
  }
}
key_freqs22 <- key_freqs22[11:20]
seasongar2fit <- lm(create_formula(key_freqs22))
y2_gar_res_new <- y2_gar_res - seasongar2fit$fitted.values
