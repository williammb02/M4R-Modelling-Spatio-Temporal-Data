# create a function that returns the index in list of coords
city_coord_index <- function(latc, lonc){
  lon_index <- which.min(abs(lon - lonc))
  lat_index <- which.min(abs(lat - latc))
  c(lon_index, lat_index)
}

# create a function that returns the index in time 
time_index <- function(t){
  match(t, climate_time)
}

# florida is bound by 31, -88, 25, -80
florida_lats = seq(from = 25, to = 31, by = 0.1)
florida_lons = seq(from = -88, to = -80, by = 0.1)

# go +- 0.8 degrees from the central coordinate

# y axis for map
miami_lats = seq(25, 26.6, 0.1)
# x axis for map
miami_lons = seq(-81, -79.4, 0.1)

tampa_lats = seq(27.2, 28.8, 0.1)
tampa_lons = seq(-83.3, -81.7, 0.1)

# create function that allows us to extract data at specific time
spatial_grid <- function(var, lats, lons, time){
  # create a matrix of 0s to store our data
  x <- matrix(0, length(lats), length(lons))
  # create the time index
  t_index <- time_index(time)
  
  for (i in 1:length(lats)){
    for (j in 1:length(lons)){
      # obtain spatial index
      spatial_index <- city_coord_index(lats[i], lons[j])
      # generate observations at each time
      x[i,j] <- import_climate_data(var, spatial_index[1], spatial_index[2])[t_index]
    }
  }
  x
}

# generate precipitation at different years
miami_tps_2020 <- spatial_grid("tp", miami_lats, miami_lons, "2020-01-01")
miami_tps_2000 <- spatial_grid("tp", miami_lats, miami_lons, "2000-01-01")
miami_tps_1980 <- spatial_grid("tp", miami_lats, miami_lons, "1980-01-01")
miami_tps_1960 <- spatial_grid("tp", miami_lats, miami_lons, "1960-01-01")

tampa_tps_1960 <- spatial_grid("tp", tampa_lats, tampa_lons, "1960-01-01")

par(mfrow = c(2,2))
plot(as.vector(miami_tps_1960))
plot(as.vector(miami_tps_1980))
plot(as.vector(miami_tps_2000))
plot(as.vector(miami_tps_2020))

# create a function that modifies data structure and fits lisa model to it for
# given latitude, longitude and variable 
lisa_vals <- function(lats, lons, var){
  x <- c()
  y <- c()
  z <- c()
  for (i in 1:length(lats)){
    for (j in 1:length(lons)){
      if (!is.na(var[i,j])) {
        x <- c(x, lats[i])
        y <- c(y, lons[j])
        z <- c(z, var[i,j])
      }
    }
  }
  cbind(x, y, z)
}
# fit lisa models to total precipitation in miami
miami_tps_1960_fit <- lisa(lisa_vals(miami_lats, miami_lons, miami_tps_1960)[,1],
                           lisa_vals(miami_lats, miami_lons, miami_tps_1960)[,2],
                           lisa_vals(miami_lats, miami_lons, miami_tps_1960)[,3], 
                           neigh = 3, resamp = 500, latlon = TRUE)
miami_tps_1980_fit <- lisa(lisa_vals(miami_lats, miami_lons, miami_tps_1980)[,1],
                           lisa_vals(miami_lats, miami_lons, miami_tps_1980)[,2],
                           lisa_vals(miami_lats, miami_lons, miami_tps_1980)[,3],
                           neigh = 3, resamp = 500, latlon = TRUE)
miami_tps_2000_fit <- lisa(lisa_vals(miami_lats, miami_lons, miami_tps_2000)[,1],
                           lisa_vals(miami_lats, miami_lons, miami_tps_2000)[,2],
                           lisa_vals(miami_lats, miami_lons, miami_tps_2000)[,3],
                           neigh = 3, resamp = 500, latlon = TRUE)
miami_tps_2020_fit <- lisa(lisa_vals(miami_lats, miami_lons, miami_tps_2020)[,1],
                           lisa_vals(miami_lats, miami_lons, miami_tps_2020)[,2],
                           lisa_vals(miami_lats, miami_lons, miami_tps_2020)[,3], 
                           neigh = 3, resamp = 500, latlon = TRUE)
# create lisa plots
par(mfrow = c(2,2))
hist(miami_tps_1960_fit$correlation, main="Miami 1960", xlab="Correlation, Moran's I")
hist(miami_tps_1980_fit$correlation, main="Miami 1980", xlab="Correlation, Moran's I")
hist(miami_tps_2000_fit$correlation, main="Miami 2000", xlab="Correlation, Moran's I")
hist(miami_tps_2020_fit$correlation, main="Miami 2020", xlab="Correlation, Moran's I")

# above mean values: red circles, below mean values: black squares
# red circles signify spatial hotspots, black squares signify spatial coldspots
par(mfrow = c(2,2))
plot(miami_tps_1960_fit, neigh.mean = FALSE, inches = 0.1, main = "LISA Plot for Miami 1960", xlab = "Latitude", ylab = "Longitude")

plot(miami_tps_1980_fit, neigh.mean = FALSE, inches = 0.1, main = "LISA Plot for Miami 1980", xlab = "Latitude", ylab = "Longitude")
plot(miami_tps_2000_fit, neigh.mean = FALSE, inches = 0.1, main = "LISA Plot for Miami 2000", xlab = "Latitude", ylab = "Longitude")
plot(miami_tps_2020_fit, neigh.mean = FALSE, inches = 0.1, main = "LISA Plot for Miami 2020", xlab = "Latitude", ylab = "Longitude")

# try a LISA plot for florida
florida_lats <- seq(25, 31, 0.1)
florida_lons <- seq(-88, -80, 0.1)

# work out coord index for the endpoints
c(city_coord_index(25, -88), city_coord_index(31, -80))

# fit LISA models to different months, and then create lisa plots
florida_tp_jul2022 <- ncvar_get(climate_data, "tp")[371:451, 191:251, 1, time_index("2022-07-01")]
x <- lisa_vals(florida_lats, florida_lons, t(florida_tp_jul2022))[,1]
y <- lisa_vals(florida_lats, florida_lons, t(florida_tp_jul2022))[,2]
z <- lisa_vals(florida_lats, florida_lons, t(florida_tp_jul2022))[,3]
florida_tp_jul2022_fit <- lisa(x, y, z, neigh = 3, resamp = 50, latlon = TRUE)

florida_tp_aug2022 <- ncvar_get(climate_data, "tp")[371:451, 191:251, 1, time_index("2022-08-01")]
x <- lisa_vals(florida_lats, florida_lons, t(florida_tp_aug2022))[,1]
y <- lisa_vals(florida_lats, florida_lons, t(florida_tp_aug2022))[,2]
z <- lisa_vals(florida_lats, florida_lons, t(florida_tp_aug2022))[,3]
florida_tp_aug2022_fit <- lisa(x, y, z, neigh = 3, resamp = 50, latlon = TRUE)

florida_tp_sep2022 <- ncvar_get(climate_data, "tp")[371:451, 191:251, 1, time_index("2022-09-01")]
x <- lisa_vals(florida_lats, florida_lons, t(florida_tp_sep2022))[,1]
y <- lisa_vals(florida_lats, florida_lons, t(florida_tp_sep2022))[,2]
z <- lisa_vals(florida_lats, florida_lons, t(florida_tp_sep2022))[,3]
florida_tp_sep2022_fit <- lisa(x, y, z, neigh = 3, resamp = 50, latlon = TRUE)

florida_tp_oct2022 <- ncvar_get(climate_data, "tp")[371:451, 191:251, 1, time_index("2022-10-01")]
x <- lisa_vals(florida_lats, florida_lons, t(florida_tp_oct2022))[,1]
y <- lisa_vals(florida_lats, florida_lons, t(florida_tp_oct2022))[,2]
z <- lisa_vals(florida_lats, florida_lons, t(florida_tp_oct2022))[,3]
florida_tp_oct2022_fit <- lisa(x, y, z, neigh = 3, resamp = 50, latlon = TRUE)

# create lisa plots for florida
plot(florida_tp_jul2022_fit, neigh.mean=FALSE, main="LISA Plot for Florida July 2022", xlab = "Latitude", ylab = "Longitude", inches=0.1)
plot(florida_tp_aug2022_fit, neigh.mean=FALSE, main="LISA Plot for Florida August 2022", xlab = "Latitude", ylab = "Longitude", inches=0.1)
plot(florida_tp_sep2022_fit, neigh.mean=FALSE, main="LISA Plot for Florida September 2022", xlab = "Latitude", ylab = "Longitude", inches=0.1)
plot(florida_tp_oct2022_fit, neigh.mean=FALSE, main="LISA Plot for Florida October 2022", xlab = "Latitude", ylab = "Longitude", inches=0.1)

us_tps_1950 <- ncvar_get(climate_data, "tp")[1:581, 1:251, 1, 1]
us_tps_2020 <- ncvar_get(climate_data, "tp")[1:581, 1:251, 1, time_index("2020-01-01")]
# note that latitude is flipped so we need to invert to plot

# create limits that can be used for all image plots for total precipitation
tpsmat <- list(us_tps_1950, us_tps_2020)
lims2 <- range(unlist(tpsmat), na.rm=TRUE)

# import skin temperature at different times
us_skt_jan_1950 <- ncvar_get(climate_data, "skt")[1:581, 1:251, 1, 1] - 273.15
us_skt_jan_2020 <- ncvar_get(climate_data, "skt")[1:581, 1:251, 1, time_index("2020-01-01")] - 273.15

us_skt_jul_1950 <- ncvar_get(climate_data, "skt")[1:581, 1:251, 1, time_index("1950-07-01")] - 273.15
us_skt_jul_2020 <- ncvar_get(climate_data, "skt")[1:581, 1:251, 1, time_index("2020-07-01")] - 273.15

# create limits that can be used for every single image plot for skin temperature 
sktmat <- list(us_skt_jan_1950, us_skt_jan_2020, us_skt_jul_1950, us_skt_jul_2020)
lims <- range(unlist(sktmat), na.rm = TRUE)

# create image plots for different variables
image.plot(lon, rev(lat), mirror.matrix(us_tps_1950), col=hcl.colors(12, "Purp"),
           xlab = "Longitude", ylab = "Latitude", main = "Total Precipitation (m) 1/1/1950",
           zlim=lims2)
image.plot(lon, rev(lat), mirror.matrix(us_tps_2020), col=hcl.colors(12, "Purp"),
           xlab = "Longitude", ylab = "Latitude", main = "Total Precipitation (m) 1/1/2020",
           zlim=lims2)
image.plot(lon, rev(lat), mirror.matrix(sktmat[[1]]), col=hcl.colors(20, "Blue-Red"),
           xlab="Longitude", ylab="Latitude", main="Skin Temperature (C) 1/1/1950",
           zlim=lims)
image.plot(lon, rev(lat), mirror.matrix(sktmat[[2]]), col=hcl.colors(20, "Blue-Red"),
           xlab="Longitude", ylab="Latitude", main="Skin Temperature (C) 1/1/2020",
           zlim=lims)
image.plot(lon, rev(lat), mirror.matrix(sktmat[[3]]), col=hcl.colors(20, "Blue-Red"),
           xlab="Longitude", ylab="Latitude", main="Skin Temperature (C) 1/7/1950",
           zlim=lims)
image.plot(lon, rev(lat), mirror.matrix(sktmat[[4]]), col=hcl.colors(20, "Blue-Red"),
           xlab="Longitude", ylab="Latitude", main="Skin Temperature (C) 1/7/2020",
           zlim=lims)
image.plot(lon, rev(lat), mirror.matrix(us_skt_jan_2020 - us_skt_jan_1950),
           col=hcl.colors(20, "Blue-Red"), xlab="Longitude", ylab="Latitude",
           main="Difference in January Skin Temperature")
image.plot(lon, rev(lat), mirror.matrix(us_skt_jul_2020 - us_skt_jul_1950),
           col=hcl.colors(20, "Blue-Red"), xlab="Longitude", ylab="Latitude",
           main="Difference in July Skin Temperature")
