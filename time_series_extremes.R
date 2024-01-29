# test for tail dependence
taildep.test(miami_w_2023, tampa_w_2023, cthresh = -0.32)
taildep.test(y_model$residuals, y2_model$residuals, cthresh = -0.41)
taildep.test(y_gar_res, y2_gar_res, cthresh = -0.4)
taildep.test(y_var_res, y2_var_res, cthresh = -0.4)

taildep.test(y_model$residuals, y2_model$residuals, cthresh = -0.05)
taildep.test(y_gar_res, y2_gar_res, cthresh = -0.05)
taildep.test(y_var_res, y2_var_res, cthresh = -0.05)

# TAIL DEPENDENCE
tde_lag <- function(u, tau){
  q1 <- quantile(miami_w_2023, u)
  q2 <- quantile(tampa_w_2023, u)
  num <- 0
  den <- 0
  l <- length(miami_w_2023)
  for (i in 1:(l-tau)){
    if ((miami_w_2023[i] > u) & (tampa_w_2023[i+tau] > u)){
      num <- num + 1
    }
    if (miami_w_2023[i] > u){
      den <- den + 1
    }
  }
  return(c(num/l, den/l, num/den))
}

tde_lag(0.8, 0)
tde_lag(0.9, 0)
tde_lag(0.95, 0)

tde_lag(0.95, 12)
tde_lag(0.95, 24)
tde_lag(0.95, 168)

tde_lag_2 <- function(x, y, u, tau){
  q1 <- quantile(x, u)
  q2 <- quantile(y, u)
  num <- 0
  den <- 0
  l <- length(x)
  for (i in 1:(l-tau)){
    if ((x[i] > u) & (y[i+tau] > u)){
      num <- num + 1
    }
    if (x[i] > u){
      den <- den + 1
    }
  }
  return(c(num/l, den/l, num/den))
}

tde_lag_2(y_model$residuals, y2_model$residuals, 0.8, 0)
tde_lag_2(y_model$residuals, y2_model$residuals, 0.9, 0)
tde_lag_2(y_model$residuals, y2_model$residuals, 0.95, 0)

tde_lag_2(y_gar_res, y2_gar_res, 0.8, 0)
tde_lag_2(y_gar_res, y2_gar_res, 0.9, 0)
tde_lag_2(y_gar_res, y2_gar_res, 0.95, 0)

# garch probabilities for all models
tde_lag_2(y_gar_res, y2_gar_res, 0.8, 0)[3]
tde_lag_2(y_gar_res, y2_gar_res, 0.9, 0)[3]
tde_lag_2(y_gar_res, y2_gar_res, 0.95, 0)[3]

tde_lag_2(y_gar_res, y3_gar_res, 0.8, 0)[3]
tde_lag_2(y_gar_res, y3_gar_res, 0.9, 0)[3]
tde_lag_2(y_gar_res, y3_gar_res, 0.95, 0)[3]

tde_lag_2(y2_gar_res, y3_gar_res, 0.8, 0)[3]
tde_lag_2(y2_gar_res, y3_gar_res, 0.9, 0)[3]
tde_lag_2(y2_gar_res, y3_gar_res, 0.95, 0)[3]
# add a lag
tde_lag_2(y_gar_res, y2_gar_res, 0.9, 24)[3]
tde_lag_2(y2_gar_res, y3_gar_res, 0.9, 24)[3]
tde_lag_2(y_gar_res, y3_gar_res, 0.9, 48)[3]

tde_lag_2(y_var_res, y2_var_res, 0.8, 0)
tde_lag_2(y_var_res, y2_var_res, 0.9, 0)
tde_lag_2(y_var_res, y2_var_res, 0.95, 0)


# CDF of ghyp distributions
q <- c(0.8, 0.9, 0.95)

set.seed(5)

pghyp(c(quantile(y_model$residuals, 0.8), quantile(y2_model$residuals, 0.8)), object = res_ghyp$best.model, lower.tail = FALSE)/0.2
pghyp(c(quantile(y_model$residuals, 0.9), quantile(y2_model$residuals, 0.9)), object = res_ghyp$best.model, lower.tail = FALSE)/0.1
pghyp(c(quantile(y_model$residuals, 0.95), quantile(y2_model$residuals, 0.95)), object = res_ghyp$best.model, lower.tail = FALSE)/0.05

pghyp(c(quantile(y_gar_res, 0.8), quantile(y2_gar_res, 0.8)), object = garch_res_ghyp$best.model, lower.tail = FALSE)/0.2
pghyp(c(quantile(y_gar_res, 0.9), quantile(y2_gar_res, 0.9)), object = garch_res_ghyp$best.model, lower.tail = FALSE)/0.1
pghyp(c(quantile(y_gar_res, 0.95), quantile(y2_gar_res, 0.95)), object = garch_res_ghyp$best.model, lower.tail = FALSE)/0.05

pghyp(c(quantile(y_var_res, 0.8), quantile(y2_var_res, 0.8)), object = var_res_ghyp$best.model, lower.tail = FALSE)/0.2
pghyp(c(quantile(y_var_res, 0.9), quantile(y2_var_res, 0.9)), object = var_res_ghyp$best.model, lower.tail = FALSE)/0.1
pghyp(c(quantile(y_var_res, 0.95), quantile(y2_var_res, 0.95)), object = var_res_ghyp$best.model, lower.tail = FALSE)/0.05








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


# return levels for GEV
rl1 <- function(p){
  return(w1e[1] - (w1e[2]/w1e[3])*(1-(-log(1-p))^(-w1e[3])))
}

rl2 <- function(p){
  return(w2e[1] - (w2e[2]/w2e[3])*(1-(-log(1-p))^(-w2e[3])))
}

# return levels for POT
zeta <- function(x, u){
  count <- 0
  n <- length(x)
  for (i in 1:n){
    if (x[i] > u){
      count <- count + 1
    }
  }
  return(count/n)
}

rlpot <- function(x, m){
  u <- quantile(x, 0.95)
  zetahat <- zeta(x, u)
  mle <- gpd.fit(x, u, npy = 24, show=FALSE)$mle
  return(u + (mle[1]/mle[2])*((m*zetahat)^(mle[2])-1))
}

c(rlpot(y_model$residuals, 24), rlpot(y_model$residuals, 168), rlpot(y2_model$residuals, 24), rlpot(y2_model$residuals, 168))
c(rlpot(y_gar_res, 24), rlpot(y_gar_res, 168), rlpot(y2_gar_res, 24), rlpot(y2_gar_res, 168))
c(rlpot(y_var_res, 24), rlpot(y_var_res, 168), rlpot(y2_var_res, 24), rlpot(y2_var_res, 168))


# extremogram 
extremogram2(cbind(y_gar_res, y2_gar_res), 0.95, 0.95, 300, 1)
extremogram2(cbind(y2_gar_res, y_gar_res), 0.95, 0.95, 300, 1)

extremogram2(cbind(y_gar_res, y3_gar_res), 0.95, 0.95, 300, 1)
extremogram2(cbind(y3_gar_res, y_gar_res), 0.95, 0.95, 300, 1)

extremogram2(cbind(y2_gar_res, y3_gar_res), 0.95, 0.95, 300, 1)
extremogram2(cbind(y3_gar_res, y2_gar_res), 0.95, 0.95, 300, 1)
