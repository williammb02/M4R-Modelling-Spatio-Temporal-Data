library(VineCopula)
library(network)
library(kdecopula)
library(rvinecopulib)
library(svines)
library(ggraph)

# method using the probability integral transform
g1 <- stepAIC.ghyp(y_gar_res, silent=TRUE)
g2 <- stepAIC.ghyp(y2_gar_res, silent=TRUE)
g3 <- stepAIC.ghyp(y3_gar_res, silent=TRUE)
g4 <- stepAIC.ghyp(y4_gar_res, silent=TRUE)
g5 <- stepAIC.ghyp(y5_gar_res, silent=TRUE)
g6 <- stepAIC.ghyp(y6_gar_res, silent=TRUE)

e1 <- pghyp(y_gar_res, object=g1$best.model)
e2 <- pghyp(y2_gar_res, object=g2$best.model)
e3 <- pghyp(y3_gar_res, object=g3$best.model)
e4 <- pghyp(y4_gar_res, object=g4$best.model)
e5 <- pghyp(y5_gar_res, object=g5$best.model)
e6 <- pghyp(y6_gar_res, object=g6$best.model)

copres2 <- cbind(e1,e2,e3,e4,e5,e6)
colnames(copres2) <- c("Miami", "Tampa", "Tallahassee", "Jacksonville", "Orlando", "Fort Myers")
cop_fit3 <- RVineStructureSelect(copres2)
cop_fit3c <- RVineStructureSelect(copres2, type = "CVine")

cop_pdf <- RVinePDF(copres2, cop_fit3)

# plots
contour(cop_fit3)
plot(cop_fit3, var_names="use")

contour(cop_fit3c)
plot(cop_fit3c, var_names="legend")

# hypothesis test of two copulas
c(cop_fit3$AIC, cop_fit3$BIC) < c(cop_fit3c$AIC, cop_fit3c$BIC)

RVineVuongTest(copres2, cop_fit3c, cop_fit3)

# simulate
set.seed(2024) 

tests <- RVineSim(5832, cop_fit3)
hist(qghyp(tests[,1], object=g1$best.model), breaks=100)


# stationary vine copula
g11 <- stepAIC.ghyp(y, silent=TRUE)
g22 <- stepAIC.ghyp(y2, silent=TRUE)
g33 <- stepAIC.ghyp(y3, silent=TRUE)
g44 <- stepAIC.ghyp(y4, silent=TRUE)
g55 <- stepAIC.ghyp(y5, silent=TRUE)
g66 <- stepAIC.ghyp(y6, silent=TRUE)

# 1 month of data
y11 <- pghyp(y[1:504], object=g11$best.model)
y22 <- pghyp(y2[1:504], object=g22$best.model)
y33 <- pghyp(y3[1:504], object=g33$best.model)
y44 <- pghyp(y4[1:504], object=g44$best.model)
y55 <- pghyp(y5[1:504], object=g55$best.model)
y66 <- pghyp(y6[1:504], object=g66$best.model)

svinedata <- cbind(y11, y22, y33, y44, y55, y66)
colnames(svinedata) <- c("M", "Tam", "Tal", "J", "O", "FM")
s_fit2 <- svinecop(svinedata, p=6)
plot(s_fit2, var_names="use")

# hourly but with precipitation too
g11t <- stepAIC.ghyp(yt, silent=TRUE)
g22t <- stepAIC.ghyp(yt2, silent=TRUE)
g33t <- stepAIC.ghyp(yt3, silent=TRUE)
g44t <- stepAIC.ghyp(yt4, silent=TRUE)
g55t <- stepAIC.ghyp(yt5, silent=TRUE)
g66t <- stepAIC.ghyp(yt6, silent=TRUE)

y11t <- pghyp(yt[1:504], object=g11t$best.model)
y11t[is.na(y11t)] <- mean(y11t[!is.na(y11t)])
y22t <- pghyp(yt2[1:504], object=g22t$best.model)
y22t[is.na(y22t)] <- mean(y22t[!is.na(y22t)])
y33t <- pghyp(yt3[1:504], object=g33t$best.model)
y33t[is.na(y33t)] <- mean(y33t[!is.na(y33t)])
y44t <- pghyp(yt4[1:504], object=g44t$best.model)
y44t[is.na(y44t)] <- mean(y44t[!is.na(y44t)])
y55t <- pghyp(yt5[1:504], object=g55t$best.model)
y55t[is.na(y55t)] <- mean(y55t[!is.na(y55t)])
y66t <- pghyp(yt6[1:504], object=g66t$best.model)
y66t[is.na(y66t)] <- mean(y66t[!is.na(y66t)])

svinedatat <- cbind(y11, y22, y33, y44, y55, y66, y11t, y22t, y33t, y44t, y55t, y66t)
colnames(svinedatat) <- c("M W", "Tam W", "Tal W", "J W", "O W", "FM W", "M P", "Tam P", "Tal P", "J P", "O P", "FM P")
s_fitt <- svinecop(svinedatat, p=3)
plot(s_fitt, var_names="use")


# hourly hurricane season, june-aug 2023



# using monthly data
g11m <- stepAIC.ghyp(ym, silent=TRUE)
g22m <- stepAIC.ghyp(ym2, silent=TRUE)
g33m <- stepAIC.ghyp(ym3, silent=TRUE)
g44m <- stepAIC.ghyp(ym4, silent=TRUE)
g55m <- stepAIC.ghyp(ym5, silent=TRUE)
g66m <- stepAIC.ghyp(ym6, silent=TRUE)

y11m <- pghyp(ym[1:480], object=g11m$best.model)
y11m[is.na(y11m)] <- mean(y11m[!is.na(y11m)])
y22m <- pghyp(ym2[1:480], object=g22m$best.model)
y22m[is.na(y22m)] <- mean(y22m[!is.na(y22m)])
y33m <- pghyp(ym3[1:480], object=g33m$best.model)
y33m[is.na(y33m)] <- mean(y33m[!is.na(y33m)])
y44m <- pghyp(ym4[1:480], object=g44m$best.model)
y44m[is.na(y44m)] <- mean(y44m[!is.na(y44m)])
y55m <- pghyp(ym5[1:480], object=g55m$best.model)
y55m[is.na(y55m)] <- mean(y55m[!is.na(y55m)])
y66m <- pghyp(ym6[1:480], object=g66m$best.model)
y66m[is.na(y66m)] <- mean(y66m[!is.na(y66m)])

svinedatam <- cbind(y11m, y22m, y33m, y44m, y55m, y66m)
colnames(svinedatam) <- c("M", "Tam", "Tal", "J", "O", "FM")
s_fitm <- svinecop(svinedatam, p=6)
plot(s_fitm, var_names="use")

# monthly but with precipitation too
g11mt <- stepAIC.ghyp(ymt, silent=TRUE)
g22mt <- stepAIC.ghyp(ymt2, silent=TRUE)
g33mt <- stepAIC.ghyp(ymt3, silent=TRUE)
g44mt <- stepAIC.ghyp(ymt4, silent=TRUE)
g55mt <- stepAIC.ghyp(ymt5, silent=TRUE)
g66mt <- stepAIC.ghyp(ymt6, silent=TRUE)

y11mt <- pghyp(ymt[1:480], object=g11m$best.model)
y11mt[is.na(y11mt)] <- mean(y11mt[!is.na(y11mt)])
y22mt <- pghyp(ymt2[1:480], object=g22m$best.model)
y22mt[is.na(y22mt)] <- mean(y22mt[!is.na(y22mt)])
y33mt <- pghyp(ymt3[1:480], object=g33m$best.model)
y33mt[is.na(y33mt)] <- mean(y33mt[!is.na(y33mt)])
y44mt <- pghyp(ymt4[1:480], object=g44m$best.model)
y44mt[is.na(y44mt)] <- mean(y44mt[!is.na(y44mt)])
y55mt <- pghyp(ymt5[1:480], object=g55m$best.model)
y55mt[is.na(y55mt)] <- mean(y55mt[!is.na(y55mt)])
y66mt <- pghyp(ymt6[1:480], object=g66m$best.model)
y66mt[is.na(y66mt)] <- mean(y66mt[!is.na(y66mt)])

svinedatamt <- cbind(y11m, y22m, y33m, y44m, y55m, y66m, y11mt, y22mt, y33mt, y44mt, y55mt, y66mt)
colnames(svinedatamt) <- c("M W", "Tam W", "Tal W", "J W", "O W", "FM W", "M P", "Tam P", "Tal P", "J P", "O P", "FM P")
s_fitmt <- svinecop(svinedatamt, p=3)
plot(s_fitmt, var_names="use")


# extreme value copula 
daily_max <- function(x){
  max_list <- split(x, rep(1:(length(x)/24), each=24))
  v <- as.vector(sapply(max_list, max))
  return(v)
}

# probability integral transform looking at extreme winds
# calculate max daily wind
w1 <- daily_max(y)
w2 <- daily_max(y2)
w3 <- daily_max(y3)
w4 <- daily_max(y4)
w5 <- daily_max(y5)
w6 <- daily_max(y6)
# fit gev dist
# loc, scale, shape
w1gev <- gev.fit(w1)
w2gev <- gev.fit(w2)
w3gev <- gev.fit(w3)
w4gev <- gev.fit(w4)
w5gev <- gev.fit(w5)
w6gev <- gev.fit(w6)

w11 <- pgev(w1, location=w1gev$mle[1], scale=w1gev$mle[2], shape=w1gev$mle[3])
w22 <- pgev(w2, location=w2gev$mle[1], scale=w2gev$mle[2], shape=w2gev$mle[3])
w33 <- pgev(w3, location=w3gev$mle[1], scale=w3gev$mle[2], shape=w3gev$mle[3])
w44 <- pgev(w4, location=w4gev$mle[1], scale=w4gev$mle[2], shape=w4gev$mle[3])
w55 <- pgev(w5, location=w5gev$mle[1], scale=w5gev$mle[2], shape=w5gev$mle[3])
w66 <- pgev(w6, location=w6gev$mle[1], scale=w6gev$mle[2], shape=w6gev$mle[3])

# calculate max daily precip
p1 <- daily_max(yt)
p2 <- daily_max(yt2)
p3 <- daily_max(yt3)
p4 <- daily_max(yt4)
p5 <- daily_max(yt5)
p6 <- daily_max(yt6)
# fit gev dist
p1gev <- gev.fit(p1)
p2gev <- gev.fit(p2)
p3gev <- gev.fit(p3)
p4gev <- gev.fit(p4)
p5gev <- gev.fit(p5)
p6gev <- gev.fit(p6)
# probability integral transform looking at extreme precipitations
p11 <- pgev(p1, location=p1gev$mle[1], scale=p1gev$mle[2], shape=p1gev$mle[3])
p22 <- pgev(p2, location=p2gev$mle[1], scale=p2gev$mle[2], shape=p2gev$mle[3])
p33 <- pgev(p3, location=p3gev$mle[1], scale=p3gev$mle[2], shape=p3gev$mle[3])
p44 <- pgev(p4, location=p4gev$mle[1], scale=p4gev$mle[2], shape=p4gev$mle[3])
p55 <- pgev(p5, location=p5gev$mle[1], scale=p5gev$mle[2], shape=p5gev$mle[3])
p66 <- pgev(p6, location=p6gev$mle[1], scale=p6gev$mle[2], shape=p6gev$mle[3])

# stationary vine copula model
svinegevdata <- cbind(w11, w22, w33, w44, w55, w66, p11, p22, p33, p44, p55, p66)
colnames(svinegevdata) <- c("M W", "Tam W", "Tal W", "J W", "O W", "FM W", "M P", "Tam P", "Tal P", "J P", "O P", "FM P")
s_fitgev <- svinecop(svinegevdata, p=3)
plot(s_fitgev, var_names="use")


# maximum daily during hurricane season, june-aug 2023
svinegevdata_h <- cbind(w11[151:243], w22[151:243], w33[151:243], w44[151:243], w55[151:243], w66[151:243], 
                        p11[151:243], p22[151:243], p33[151:243], p44[151:243], p55[151:243], p66[151:243])
colnames(svinegevdata_h) <- c("M W", "Tam W", "Tal W", "J W", "O W", "FM W", "M P", "Tam P", "Tal P", "J P", "O P", "FM P")
s_fitgev_h <- svinecop(svinegevdata_h, p=3)
plot(s_fitgev_h, var_names="use")




# plots to show the dependence 
flo <- map_data("state", region = "florida")
cities <- data.frame(
  city = c("Miami (1)", "Tampa (2)", "Tallahassee (3)", "Jacksonville (4)", "Orlando (5)", "Fort Myers (6)"),
  lon = c(-80.2, -82.5, -84.2, -81.7, -81.4, -81.9),
  lat = c(25.8, 28.0, 30.4, 30.3, 28.5, 26.6)
)
edgesd <- data.frame(
  from = c("Miami (1)", "Fort Meyers (6)", "Tampa (2)", "Orlando (5)", "Jacksonville (4)"),
  to = c("Fort Meyers (6)", "Tampa(2)", "Orlando (5)", "Jacksonville (4)", "Tallahassee (3)"), 
  lon = c(-80.2, -81.9, -82.5, -81.4, -81.7),
  lat = c(25.8, 26.6, 28.0, 28.5, 30.3),
  lon2 = c(-81.9, -82.5, -81.4, -81.7, -84.2),
  lat2 = c(26.6, 28.0, 28.5, 30.3, 30.4)
)
edgesc <- data.frame(
  from = c("Orlando (5)", "Orlando (5)", "Orlando (5)", "Orlando (5)", "Orlando (5)"),
  to = c("Miami (1)", "Tampa (2)", "Tallahassee (3)", "Jacksonville (4)", "Fort Meyers (6)"),
  lon = c(-81.4, -81.4, -81.4, -81.4, -81.4),
  lat = c(28.5, 28.5, 28.5, 28.5, 28.5),
  lon2 = c(-80.2, -82.5, -84.2, -81.7, -81.9),
  lat2 = c(25.8, 28.0, 30.4, 30.3, 26.6)
)
edgess1 <- data.frame(
  from = c("Fort Meyers", "Miami", "Jacksonville", "Jacksonville", "Tampa"),
  to = c("Miami", "Jacksonville", "Tallahassee", "Tampa", "Orlando"),
  lon = c(-81.9, -80.2, -81.7, -81.7, -82.5),
  lat = c(26.6, 25.8, 30.3, 30.3, 28.0),
  lon2 = c(-80.2, -81.7, -84.2, -82.5, -81.4),
  lat2 = c(25.8, 30.3, 30.4, 28.0, 28.5)
)
edgess2 <- data.frame(
  from = c("Miami", "Tampa", "Tampa", "Tampa", "Orlando"),
  to = c("Tampa", "Tallahassee", "Fort Myers", "Orlando", "Jacksonville"),
  lon = c(-80.2, -82.5, -82.5, -82.5, -81.4),
  lat = c(25.8, 28.0, 28.0, 28.0, 28.5),
  lon2 = c(-82.5, -84.2, -81.9, -81.4, -81.7),
  lat2 = c(28.0, 30.4, 26.6, 28.5, 30.3)
)
edgess3 <- data.frame(
  from = c("Miami", "Fort Meyers", "Tampa", "Orlando", "Jacksonville"),
  to = c("Fort Meyers", "Tampa", "Orlando", "Jacksonville", "Tallahassee"), 
  lon = c(-80.2, -81.9, -82.5, -81.4, -81.7),
  lat = c(25.8, 26.6, 28.0, 28.5, 30.3),
  lon2 = c(-81.9, -82.5, -81.4, -81.7, -84.2),
  lat2 = c(26.6, 28.0, 28.5, 30.3, 30.4)
)
edgestphour <- data.frame(
  from = c("Miami", "Orlando", "Tallahassee", "Tallahassee"),
  to = c("Orlando", "Tallahassee", "Jacksonville", "Tampa"), 
  lon = c(-80.2, -81.4, -84.2, -84.2),
  lat = c(25.8, 28.5, 30.4, 30.4),
  lon2 = c(-81.4, -84.2, -81.7, -82.5),
  lat2 = c(28.5, 30.4, 30.3, 28.0)
)
edgesextra <- data.frame(
  from = c("Tallahassee"),
  to = c("Jacksonville"),
  lon = c(-84.2),
  lat = c(30.4),
  lon2 = c(-81.7),
  lat2 = c(30.3)
)
edgestpmonth <- data.frame(
  from = c("Miami", "Jacksonville", "Tampa", "Tampa", "Tampa"),
  to = c("Jacksonville", "Tampa", "Orlando", "Fort Myers", "Tallahassee"), 
  lon = c(-80.2, -81.7, -82.5, -82.5, -82.5),
  lat = c(25.8, 30.3, 28.0, 28.0, 28.0),
  lon2 = c(-81.7, -82.5, -81.4, -81.9, -84.2),
  lat2 = c(30.3, 28.0, 28.5, 26.6, 30.4)
)
edgesextram <- data.frame(
  from = c("Tampa", "Tampa"),
  to = c("Fort Myers", "Orlando"),
  lon = c(-82.5, -82.5),
  lat = c(28.0, 28.0),
  lon2 = c(-81.9, -81.4),
  lat2 = c(26.6, 28.5)
)
edgestpgev <- data.frame(
  from = c("Miami", "Jacksonville", "Tampa", "Tallahassee"),
  to = c("Fort Meyers", "Tampa", "Tallahassee", "Orlando"), 
  lon = c(-80.2, -81.7, -82.5, -84.2),
  lat = c(25.8, 30.3, 28.0, 30.4),
  lon2 = c(-81.9, -82.5, -84.2, -81.4),
  lat2 = c(26.6, 28.0, 30.4, 28.5)
)
edgesextragev<- data.frame(
  from = c("Miami"),
  to = c("Fort Myers"),
  lon = c(-80.2),
  lat = c(25.8),
  lon2 = c(-81.9),
  lat2 = c(26.6)
)
linkh <- data.frame(
  from = c("Orlando", "Tallahassee"),
  to = c("Fort Myers", "Miami"),
  lon = c(-81.4, -84.2),
  lat = c(28.5, 30.4),
  lon2 = c(-81.9, -80.2),
  lat2 = c(26.6, 25.8)
)
linkm <- data.frame(
  from = c("Tallahassee"),
  to = c("Fort Myers"),
  lon = c(-84.2),
  lat = c(30.4),
  lon2 = c(-81.9),
  lat2 = c(26.6)
)
linkgev <- data.frame(
  from = c("Tampa", "Miami"),
  to = c("Jacksonville", "Miami"),
  lon = c(-82.5, -80.2),
  lat = c(28.0, 25.8),
  lon2 = c(-81.7, -80.2),
  lat2 = c(30.3, 25.8)
)

edgeswhurricane <- data.frame(
  from = c("Miami", "Fort Myers", "Fort Myers", "Orlando", "Jacksonville"),
  to = c("Fort Myers", "Tampa", "Orlando", "Jacksonville", "Tallahassee"), 
  lon = c(-80.2, -81.9, -81.9, -81.4, -81.7),
  lat = c(25.8, 26.6, 26.6, 28.5, 30.3),
  lon2 = c(-81.9, -82.5, -81.4, -81.7, -84.2),
  lat2 = c(26.6, 28.0, 28.5, 30.3, 30.4)
)
edgestphurricane <- data.frame(
  from = c("Jacksonville", "Tampa", "Tampa", "Tampa"),
  to = c("Tampa", "Orlando", "Fort Myers", "Tallahassee"), 
  lon = c(-81.7, -82.5, -82.5, -82.5),
  lat = c(30.3, 28.0, 28.0, 28.0),
  lon2 = c(-82.5, -81.4, -81.9, -84.2),
  lat2 = c(28.0, 28.5, 26.6, 30.4)
)

ggplot() +
  ggtitle("Map of Florida and Cities in the D-Vine and S-Vine Copula Models") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_segment(data = edgesd, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 1) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)
  
ggplot() +
  ggtitle("Map of Florida and Cities in the C-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_segment(data = edgesc, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 1) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)

# hourly joint plot
ggplot() +
  ggtitle("Map of Florida and Cities in the S-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "black", size = 3) +
  geom_segment(data = edgesd, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 0.5) +
  geom_segment(data = edgestphour, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "blue", size = 0.5) +
  geom_segment(data = linkh, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "black", size = 0.5, linetype=2) +
  geom_segment(data = edgesextra,
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "purple", size = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)
# monthly joint plot
ggplot() +
  ggtitle("Map of Florida and Cities in the S-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "black", size = 3) +
  geom_segment(data = edgestpmonth, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "blue", size = 0.5) +
  geom_segment(data = edgesd, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 0.5) +
  geom_segment(data = linkm, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "black", size = 0.5, linetype=2) +
  geom_segment(data = edgesextram, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "purple", size = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)
# daily maximum joint plot 
ggplot() +
  ggtitle("Map of Florida and Cities in the S-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "black", size = 3) +
  geom_segment(data = edgestpgev, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "blue", size = 0.5) +
  geom_segment(data = edgesd, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 0.5) +
  geom_segment(data = linkgev, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "black", size = 0.8, linetype=2) +
  geom_segment(data = edgesextragev, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "purple", size = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)

# hurricane season wind, then precip
ggplot() +
  ggtitle("Map of Florida and Cities in the S-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_segment(data = edgeswhurricane, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "red", size = 1) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)

ggplot() +
  ggtitle("Map of Florida and Cities in the S-Vine Copula Model") +
  geom_polygon(data = flo, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "black") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "blue", size = 3) +
  geom_segment(data = edgestphurricane, 
               aes(x = lon, y = lat, xend = lon2, yend = lat2, group = NULL),
               color = "blue", size = 1) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), 
            color = "black", size = 4, vjust = -1) +
  coord_fixed(ratio = 1.3)


