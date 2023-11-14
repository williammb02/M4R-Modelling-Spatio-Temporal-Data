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


# tail dependence coefficients
taildeps <- function(x, y){
  lambda <- c()
  t <- seq(0.05, 0.95, by=0.025)
  for(i in 1:length(t)){
    lambda[i] <- taildep(x, y, t[i], type="chi")
  }
  lambda
}

t <- seq(0.05, 0.95, by=0.025)

# components of wind speed at same location
plot(t, taildeps(miami_u10, miami_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(tampa_u10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(tallahassee_u10, tallahassee_v10), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# same component, different location
plot(t, taildeps(miami_u10, tampa_u10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa u component")
plot(t, taildeps(abs(miami_u10), abs(tampa_u10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(miami_v10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa v component")
plot(t, taildeps(abs(miami_v10), abs(tampa_v10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
# ones below are particularly notable, consider absolute
plot(t, taildeps(tallahassee_u10, tampa_u10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa u component")
plot(t, taildeps(abs(tallahassee_u10), abs(tampa_u10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))
plot(t, taildeps(tallahassee_v10, tampa_v10), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa v component")
plot(t, taildeps(abs(tallahassee_v10), abs(tampa_v10)), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# precipitation
plot(t, taildeps(tallahassee_tp, tampa_tp), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))

# wind speed
plot(t, taildeps(miami_speed, tampa_speed), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Miami and Tampa speed")
plot(t, taildeps(talla_speed, tampa_speed), xlab="Probability Treshold", ylab="Coefficient", 
     ylim=c(0,1), main="Tallahassee and Tampa speed")
plot(t, taildeps(talla_speed, miami_speed), xlab="Probability Treshold", ylab="Coefficient", ylim=c(0,1))


# TAIL DEPENDENCE WITH NEGATIVE VERSION OF U
# TAIL DEPENDENCE WITH REMAINDERS OF TIME SERIES

# extremogram analysis
# we need stationary time series to use the extremogram
# use remainders from temporal.R

extremogram2(cbind(miami_u10_rem, miami_v10_rem), 0.95, 0.95, 100, 1)
extremogram2(cbind(miami_u10_rem, tampa_u10_rem), 0.95, 0.95, 100, 1)
extremogram2(cbind(miami_u10_rem, miami_tp_rem), 0.95, 0.95, 100, 1)




uvm <- cbind(miami_u10, miami_v10)
utm <- cbind(miami_u10, miami_tp)
t_mt <- cbind(miami_tp, tampa_tp)
mtp_s <- cbind(miami_speed, tampa_speed)
mth_s <- cbind(miami_speed, talla_speed)
tt_s <- cbind(tampa_speed, talla_speed)
# nice extremal dependence
extremogram2(uvm, 0.95, 0.05, 100, 4)
extremogram2(utm, 0.95, 0.95, 100, 1)
extremogram2(utm, 0.95, 0.05, 100, 4)
extremogram2(t_mt, 0.95, 0.95, 100, 1)
extremogram2(t_mt, 0.05, 0.05, 100, 2)
extremogram2(t_mt, 0.05, 0.95, 100, 3)
extremogram2(t_mt, 0.95, 0.05, 100, 4)

extremogram1(miami_u10, 0.95, 100, 1)
extremogram1(miami_u10, 0.05, 100, 2)
# stronger pattern in v
extremogram1(miami_v10, 0.95, 100, 1)
extremogram1(miami_v10, 0.05, 100, 2)

# one below is very interesting
extremogram2(tt_s, 0.95, 0.95, 100, 1)
extremogram2(tt_s, 0.05, 0.05, 100, 2)

