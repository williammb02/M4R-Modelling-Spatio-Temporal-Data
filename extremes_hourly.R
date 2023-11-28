# hourly data analysis
hourly_miami_tp <- stl_rem_hourly(miami_tp_2023)
hourly_miami_u <- stl_rem_hourly(miami_u10_2023)
hourly_miami_v <- stl_rem_hourly(miami_v10_2023)

hourly_tampa_tp <- stl_rem_hourly(tampa_tp_2023)
hourly_tampa_u <- stl_rem_hourly(tampa_u10_2023)
hourly_tampa_v <- stl_rem_hourly(tampa_v10_2023)

hourly_talla_tp <- stl_rem_hourly(talla_tp_2023)
hourly_talla_u <- stl_rem_hourly(talla_u10_2023)
hourly_talla_v <- stl_rem_hourly(talla_v10_2023)

miami_w_2023 <- sqrt(miami_u10_2023^2 + miami_v10_2023^2)
tampa_w_2023 <- sqrt(tampa_u10_2023^2 + tampa_v10_2023^2)
talla_w_2023 <- sqrt(talla_u10_2023^2 + talla_v10_2023^2)

hourly_miami_w <- stl_rem_hourly(miami_w_2023)
hourly_tampa_w <- stl_rem_hourly(tampa_w_2023)
hourly_talla_w <- stl_rem_hourly(talla_w_2023)

# tail dependence coefficients
taildeps <- function(t, x, y){
  lambda <- c()
  for(i in 1:length(t)){
    lambda[i] <- taildep(x, y, t[i], type="chi")
  }
  lambda
}

t <- seq(0.05, 0.95, by=0.025)
t2 <- seq(0.7, 0.99, by = 0.01)

# tail dependence coefficients
# u and v components
plot(t2, taildeps(t2, hourly_miami_u, hourly_tampa_u), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa u component remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_u, hourly_tampa_u), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa u component remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_u, hourly_miami_u), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Miami u component remainder", type="l")
plot(t2, taildeps(t2, hourly_miami_v, hourly_tampa_v), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa v component remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_v, hourly_tampa_v), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa v component remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_v, hourly_miami_v), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Miami v component remainder", type="l")
# talla and miami not much dependence
# miami and tampa / tampa and talla dependence, they are closer

# total precipitation
plot(t2, taildeps(t2, hourly_miami_tp, hourly_tampa_tp), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa precipitation remainder", type="l")
plot(t2, taildeps(t2, hourly_miami_tp, hourly_talla_tp), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tallahassee precipitation remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_tp, hourly_tampa_tp), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa precipitation remainder", type="l")
# talla and tampa have some dependence, others not so much 


# wind speed
plot(t2, taildeps(t2, hourly_tampa_w, hourly_miami_w), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa wind speed remainder", type="l")
plot(t2, taildeps(t2, hourly_talla_w, hourly_miami_w), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tallahassee wind speed remainder", type="l")
plot(t2, taildeps(t2, hourly_tampa_w, hourly_talla_w), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa wind speed remainder", type="l")
# strong dependence in the limit, good source of extremal analysis

plot(t2, taildeps(t2, miami_u10_2023, tampa_u10_2023), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa u component", type="l")
plot(t2, taildeps(t2, miami_tp_2023, tampa_tp_2023), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa precipitation", type="l")
plot(t2, taildeps(t2, miami_w_2023, tampa_w_2023), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa wind speed", type="l")


# testing tail dependence against tail independence
# look for between 10-15% exceedance rate
# p value low means independent, p value high means dependent

# wind components
taildep.test(miami_u10_2023, miami_v10_2023, cthresh = -0.50)
taildep.test(miami_u10_2023, miami_v10_2023, cthresh = -0.42)

taildep.test(miami_u10_2023, tampa_u10_2023, cthresh = -0.40)
taildep.test(miami_u10_2023, tampa_u10_2023, cthresh = -0.32)

taildep.test(miami_v10_2023, tampa_v10_2023, cthresh = -0.38)
taildep.test(miami_v10_2023, tampa_v10_2023, cthresh = -0.26)

# wind speed
taildep.test(miami_w_2023, tampa_w_2023, cthresh = -0.34)
taildep.test(miami_w_2023, talla_w_2023, cthresh = -0.37)
taildep.test(tampa_w_2023, talla_w_2023, cthresh = -0.28)

# remainders
taildep.test(hourly_miami_v, hourly_tampa_v, cthresh = -0.27)
taildep.test(hourly_miami_w, hourly_tampa_w, cthresh = -0.33)


# extremogram analysis
# more insight that monthly
extremogram2(cbind(hourly_miami_tp, hourly_tampa_tp), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_tampa_tp, hourly_miami_tp), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_miami_v, hourly_miami_u), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_miami_u, hourly_miami_v), 0.95, 0.95, 300, 1)

# wind speed - more dependence within the next day or two (up to lag 50 hours)
# again varies spatially 
extremogram2(cbind(hourly_miami_w, hourly_tampa_w), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_talla_w, hourly_tampa_w), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_talla_w, hourly_miami_w), 0.95, 0.95, 300, 1)

extremogram2(cbind(hourly_tampa_w, hourly_miami_w), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_tampa_w, hourly_talla_w), 0.95, 0.95, 300, 1)
extremogram2(cbind(hourly_miami_w, hourly_talla_w), 0.95, 0.95, 300, 1)



extremogram2(cbind(hourly_miami_tp, hourly_tampa_tp), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_tampa_tp, hourly_miami_tp), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_miami_v, hourly_miami_u), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_miami_u, hourly_miami_v), 0.9, 0.9, 300, 1)

extremogram2(cbind(hourly_miami_w, hourly_tampa_w), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_talla_w, hourly_tampa_w), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_talla_w, hourly_miami_w), 0.9, 0.9, 300, 1)

extremogram2(cbind(hourly_tampa_w, hourly_miami_w), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_tampa_w, hourly_talla_w), 0.9, 0.9, 300, 1)
extremogram2(cbind(hourly_miami_w, hourly_talla_w), 0.9, 0.9, 300, 1)

