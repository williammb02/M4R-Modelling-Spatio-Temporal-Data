# test for tail dependence
taildep.test(miami_w_2023, tampa_w_2023, cthresh = -0.32)
taildep.test(y_model$residuals, y2_model$residuals, cthresh = -0.41)
taildep.test(y_gar_res, y2_gar_res, cthresh = -0.4)
taildep.test(y_var_res, y2_var_res, cthresh = -0.4)

taildep.test(y_model$residuals, y2_model$residuals, cthresh = -0.05)
taildep.test(y_gar_res, y2_gar_res, cthresh = -0.05)
taildep.test(y_var_res, y2_var_res, cthresh = -0.05)


# look at maximum daily rainfall
daily_max <- function(x){
  max_list <- split(x, rep(1:(length(x)/24), each=24))
  v <- as.vector(sapply(max_list, max))
  return(v)
}

w1 <- daily_max(miami_w_2023)
w2 <- daily_max(tampa_w_2023)

w1gev <- gev.fit(w1)
w2gev <- gev.fit(w2)
w1e <- w1gev$mle
w2e <- w2gev$mle

w1gpd <- gpd.fit(miami_w_2023, quantile(miami_w_2023, 0.95), npy=24)
w2gpd <- gpd.fit(tampa_w_2023, quantile(tampa_w_2023, 0.95), npy=24)

# univariate cdfs
pot_dist1 <- function(x){
  return(1-(1 + w1gpd$mle[2]*(x-quantile(miami_w_2023, 0.95))/w1gpd$mle[1])^(-1/w1gpd$mle[2]))
}

pot_dist2 <- function(x){
  return(1-(1 + w2gpd$mle[2]*(x-quantile(tampa_w_2023, 0.95))/w2gpd$mle[1])^(-1/w2gpd$mle[2]))
}

# bivariate cdf
bpot_dist <- function(x, y, alpha){
  x1 <- -(log(pot_dist1(x)))^(-1)
  y1 <- -(log(pot_dist2(y)))^(-1)
  v <- (x1^(-1/alpha) + y1^(-1/alpha))^alpha
  return(exp(-v))
}

# bivariate survival function
bpot_surv <- function(x, y, alpha){
  return(1-pot_dist1(x)-pot_dist2(y)+bpot_dist(x, y, alpha))
}

# tail dependence
bpot_surv(quantile(miami_w_2023, 0.99), quantile(tampa_w_2023, 0.99), 0.5)/(1-pot_dist1(quantile(miami_w_2023, 0.99)))


# POT models for the remainders and estimates of exceedance probabilities

r11gpd <- gpd.fit(y_model$residuals, quantile(y_model$residuals, 0.95), npy=24)
r12gpd <- gpd.fit(y2_model$residuals, quantile(y2_model$residuals, 0.95), npy=24)

r21gpd <- gpd.fit(y_gar_res, quantile(y_gar_res, 0.95), npy=24)
r22gpd <- gpd.fit(y2_gar_res, quantile(y2_gar_res, 0.95), npy=24)

r31gpd <- gpd.fit(y_var_res, quantile(y_var_res, 0.95), npy=24)
r32gpd <- gpd.fit(y2_var_res, quantile(y2_var_res, 0.95), npy=24)

pot_dist_rem1 <- function(x, i){
  if(i == 1) {
    return(1-(1 + r11gpd$mle[2]*(x-quantile(y_model$residuals, 0.95))/r11gpd$mle[1])^(-1/r11gpd$mle[2]))
  } else if(i == 2){
    return(1-(1 + r21gpd$mle[2]*(x-quantile(y_gar_res, 0.95))/r21gpd$mle[1])^(-1/r21gpd$mle[2]))
  } else if(i == 3){
    return(1-(1 + r31gpd$mle[2]*(x-quantile(y_var_res, 0.95))/r31gpd$mle[1])^(-1/r31gpd$mle[2]))
  }
}

pot_dist_rem2 <- function(x, i){
  if(i == 1) {
    return(1-(1 + r12gpd$mle[2]*(x-quantile(y_model$residuals, 0.95))/r12gpd$mle[1])^(-1/r12gpd$mle[2]))
  } else if(i == 2) {
    return(1-(1 + r22gpd$mle[2]*(x-quantile(y_gar_res, 0.95))/r22gpd$mle[1])^(-1/r22gpd$mle[2]))
  } else if(i == 3) {
    return(1-(1 + r32gpd$mle[2]*(x-quantile(y_var_res, 0.95))/r32gpd$mle[1])^(-1/r32gpd$mle[2]))
  }
}

bpot_dist_rem <- function(x, y, i, alpha){
  x1 <- -(log(pot_dist_rem1(x, i)))^(-1)
  y1 <- -(log(pot_dist_rem2(y, i)))^(-1)
  v <- (x1^(-1/alpha) + y1^(-1/alpha))^alpha
  return(exp(-v))
}

bpot_surv_rem <- function(x, y, i, alpha){
  return(1-pot_dist_rem1(x, i)-pot_dist_rem2(y, i)+bpot_dist_rem(x, y, i, alpha))
}


bpot_surv_rem(quantile(y_model$residuals, 0.99), quantile(y2_model$residuals, 0.99), 1, 0.5)/(1-pot_dist_rem1(quantile(y_model$residuals, 0.99), 1))
bpot_surv_rem(quantile(y_gar_res, 0.99), quantile(y2_gar_res, 0.99), 2, 0.5)/(1-pot_dist_rem1(quantile(y_gar_res, 0.99), 2))
bpot_surv_rem(quantile(y_var_res, 0.99), quantile(y2_var_res, 0.99), 3, 0.5)/(1-pot_dist_rem1(quantile(y_var_res, 0.99), 3))


