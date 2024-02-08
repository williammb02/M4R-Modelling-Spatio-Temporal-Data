library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ismev)
library(sp)
library(sf)
library(maps)
library(mapdata)
library(SynchWave)
library(Thermimage)
library(ncf)
library(GCor)
library(ghyp)
library(forecast)
library(feasts)
library(fitdistrplus)
library(FAdist)
library(mclust)
library(extremogram)
library(mixR)
library(extRemes)
library(rugarch)
library(vars)
library(VineCopula)

# read in the data
climate_data <- nc_open('m4rdata.nc')
lon <- ncvar_get(climate_data, "longitude")
lat <- ncvar_get(climate_data, "latitude")
time <- ncvar_get(climate_data, "time")

# time is given in hours since 01/01/1990 so we need to transform the time
# this now numbers it by date.

climate_time <- nc.get.time.series(climate_data, time.dim.name = "time")

# create a function that imports data appropriately from the file by giving appropriate lat and lon coordinates
import_climate_data <- function(var, lon_index, lat_index){
  beforedata <- nc.get.var.subset.by.axes(climate_data, var, axis.indices = list(X = lon_index, Y = lat_index))
  newdata <- as.vector(beforedata)[c(TRUE, FALSE)]
  data <- newdata[1:882]
  data
}

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

# create a function that returns a temporal periodogram of a time series
my_periodogram <- function(x){
  N <- length(x)
  S <- (1/N)*abs(fft(x - mean(x)))^2
  fftshift(S)
}

# a collection if indexes for cities
nyc_index <- city_coord_index(40.7, -74.0)
la_index <- city_coord_index(34.1, -118.2)
austin_index <- city_coord_index(30.3, -97.7)
dc_index <- city_coord_index(38.9, -77.0)
miami_index <- city_coord_index(25.8, -80.2)
tampa_index <- city_coord_index(28.0, -82.5)
tallahassee_index <- city_coord_index(30.4, -84.2)


# read in new hourly data
miami_tp_2023_nc <- nc_open('miami_tp_2023.nc')
miami_u10_2023_nc <- nc_open('miami_u10_2023.nc')
miami_v10_2023_nc <- nc_open('miami_v10_2023.nc')

miami_tp_2023 <- as.vector(nc.get.var.subset.by.axes(miami_tp_2023_nc, "tp", axis.indices=list(X=1,Y=1)))
miami_u10_2023 <- as.vector(nc.get.var.subset.by.axes(miami_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
miami_v10_2023 <- as.vector(nc.get.var.subset.by.axes(miami_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

tampa_tp_2023_nc <- nc_open('tampa_tp_2023.nc')
tampa_u10_2023_nc <- nc_open('tampa_u10_2023.nc')
tampa_v10_2023_nc <- nc_open('tampa_v10_2023.nc')

tampa_tp_2023 <- as.vector(nc.get.var.subset.by.axes(tampa_tp_2023_nc, "tp", axis.indices=list(X=1,Y=1)))
tampa_u10_2023 <- as.vector(nc.get.var.subset.by.axes(tampa_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
tampa_v10_2023 <- as.vector(nc.get.var.subset.by.axes(tampa_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

talla_tp_2023_nc <- nc_open('talla_tp_2023.nc')
talla_u10_2023_nc <- nc_open('talla_u10_2023.nc')
talla_v10_2023_nc <- nc_open('talla_v10_2023.nc')

talla_tp_2023 <- as.vector(nc.get.var.subset.by.axes(talla_tp_2023_nc, "tp", axis.indices=list(X=1,Y=1)))
talla_u10_2023 <- as.vector(nc.get.var.subset.by.axes(talla_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
talla_v10_2023 <- as.vector(nc.get.var.subset.by.axes(talla_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

jack_u10_2023_nc <- nc_open('jack_u10_2023.nc')
jack_v10_2023_nc <- nc_open('jack_v10_2023.nc')
jack_u10_2023 <- as.vector(nc.get.var.subset.by.axes(jack_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
jack_v10_2023 <- as.vector(nc.get.var.subset.by.axes(jack_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

orla_u10_2023_nc <- nc_open('orla_u10_2023.nc')
orla_v10_2023_nc <- nc_open('orla_v10_2023.nc')
orla_u10_2023 <- as.vector(nc.get.var.subset.by.axes(orla_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
orla_v10_2023 <- as.vector(nc.get.var.subset.by.axes(orla_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

fort_u10_2023_nc <- nc_open('fort_u10_2023.nc')
fort_v10_2023_nc <- nc_open('fort_v10_2023.nc')
fort_u10_2023 <- as.vector(nc.get.var.subset.by.axes(fort_u10_2023_nc, "u10", axis.indices=list(X=1,Y=1)))
fort_v10_2023 <- as.vector(nc.get.var.subset.by.axes(fort_v10_2023_nc, "v10", axis.indices=list(X=1,Y=1)))

# visualise location of all the cities in Florida
flo <- map_data("state", region = "florida")
cities <- data.frame(
  city = c("Miami (1)", "Tampa (2)", "Tallahassee (3)", "Jacksonville (4)", "Orlando (5)", "Fort Myers (6)"),
  lon = c(-80.2, -82.5, -84.2, -81.7, -81.4, -81.9),
  lat = c(25.8, 28.0, 30.4, 30.3, 28.5, 26.6)
)

ggplot() +
  ggtitle("Map of Florida and Cities in the Copula Models") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3) +
  theme_void()

