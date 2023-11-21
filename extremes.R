miami_tp <- import_climate_data("tp", miami_index[1], miami_index[2])
tampa_tp <- import_climate_data("tp", tampa_index[1], tampa_index[2])
tallahassee_tp <- import_climate_data("tp", tallahassee_index[1], tallahassee_index[2])

miami_u10 <- import_climate_data("u10", miami_index[1], miami_index[2])
tampa_u10 <- import_climate_data("u10", tampa_index[1], tampa_index[2])
tallahassee_u10 <- import_climate_data("u10", tallahassee_index[1], tallahassee_index[2])

miami_v10 <- import_climate_data("v10", miami_index[1], miami_index[2])
tampa_v10 <- import_climate_data("v10", tampa_index[1], tampa_index[2])
tallahassee_v10 <- import_climate_data("v10", tallahassee_index[1], tallahassee_index[2])

miami_speed <- sqrt(miami_u10^2 + miami_v10^2)
tampa_speed <- sqrt(tampa_u10^2 + tampa_v10^2)
talla_speed <- sqrt(tallahassee_u10^2 + tallahassee_v10^2)

miami_w_rem <- stl_rem(miami_speed)
tampa_w_rem <- stl_rem(tampa_speed)
talla_w_rem <- stl_rem(talla_speed)

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

# components of wind speed at same location
plot(t, taildeps(t, miami_u10, miami_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(t, tampa_u10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(t, tallahassee_u10, tallahassee_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# same component, different location
plot(t, taildeps(t, miami_u10, tampa_u10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa u component")
plot(t, taildeps(t, abs(miami_u10), abs(tampa_u10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(t, miami_v10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa v component")
plot(t, taildeps(t, abs(miami_v10), abs(tampa_v10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
# ones below are particularly notable, consider absolute
plot(t, taildeps(t, tallahassee_u10, tampa_u10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa u component")
plot(t, taildeps(t, abs(tallahassee_u10), abs(tampa_u10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(t, tallahassee_v10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa v component")
plot(t, taildeps(t, abs(tallahassee_v10), abs(tampa_v10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# precipitation
plot(t, taildeps(t, tallahassee_tp, tampa_tp), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# wind speed
plot(t, taildeps(t, miami_speed, tampa_speed), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa speed")
plot(t, taildeps(t, talla_speed, tampa_speed), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa speed")
plot(t, taildeps(t, talla_speed, miami_speed), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))


# negative version of u for miami
# completely independent
plot(t, taildeps(-miami_u10, miami_v10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami -u and v component")

# remainders of time series
plot(t, taildeps(miami_u10_rem, tampa_u10_rem), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa u component remainder")
plot(t, taildeps(miami_v10_rem, tampa_v10_rem), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa v component remainder")
# ones below are particularly notable, consider absolute
plot(t, taildeps(talla_u10_rem, tampa_u10_rem), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa u component remainder")
plot(t, taildeps(talla_v10_rem, tampa_v10_rem), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa v component remainder")


# extremogram analysis
# we need stationary time series to use the extremogram
# use remainders from temporal.R

# most 'extreme' ones
extremogram2(cbind(tampa_v10, tampa_tp), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_v10_rem, tampa_tp_rem), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_u10, tampa_tp), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_u10_rem, tampa_tp_rem), 0.95, 0.95, 100, 1)

extremogram2(cbind(tampa_tp, tampa_u10), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_tp_rem, tampa_u10_rem), 0.95, 0.95, 100, 1)

# nice extremal dependence
extremogram2(cbind(miami_u10, miami_v10), 0.05, 0.95, 100, 3)

extremogram1(miami_u10, 0.95, 100, 1)
extremogram1(miami_u10, 0.05, 100, 2)
# stronger pattern in v, good extremal dependence
extremogram1(miami_v10, 0.95, 100, 1)
extremogram1(miami_v10, 0.05, 100, 2)


# one below is very interesting
extremogram2(cbind(tampa_speed, talla_speed), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_w_rem, talla_w_rem), 0.95, 0.95, 100, 1)

extremogram2(cbind(talla_speed, tampa_speed), 0.95, 0.95, 100, 1)
extremogram2(cbind(talla_w_rem, tampa_w_rem), 0.95, 0.95, 100, 1)
# not as extreme
extremogram2(cbind(talla_w_rem, miami_w_rem), 0.95, 0.95, 100, 1)
extremogram2(cbind(tampa_w_rem, miami_w_rem), 0.95, 0.95, 100, 1)




# hourly data analysis
hourly_miami_tp <- stl_rem_hourly(miami_tp_2023)
hourly_miami_u <- stl_rem_hourly(miami_u10_2023)
hourly_miami_v <- stl_rem_hourly(miami_v10_2023)

extremogram2(cbind(hourly_miami_tp, hourly_miami_u), 0.95, 0.95, 100, 1)
extremogram2(cbind(hourly_miami_u, hourly_miami_tp), 0.95, 0.95, 100, 1)
extremogram2(cbind(hourly_miami_v, hourly_miami_u), 0.95, 0.95, 100, 1)
extremogram2(cbind(hourly_miami_u, hourly_miami_v), 0.95, 0.95, 100, 1)




