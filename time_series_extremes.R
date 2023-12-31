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
