daily_max <- function(x){
  max_list <- split(x, rep(1:(length(x)/24), each=24))
  v <- as.vector(sapply(max_list, max))
  return(v)
}

w1 <- daily_max(miami_w_2023)
w2 <- daily_max(tampa_w_2023)

w1gev <- gev.fit(w1)
w2gev <- gev.fit(w2)

w1gpd <- gpd.fit(miami_w_2023, quantile(miami_w_2023, 0.95), npy=24)
w2gpd <- gpd.fit(tampa_w_2023, quantile(tampa_w_2023, 0.95), npy=24)
