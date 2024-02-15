library(VineCopula)
library(network)
library(kdecopula)
library(rvinecopulib)
library(svines)
library(ggraph)

# fit a copula model to residuals at 3 locations

cop_res <- cbind(y_gar_res, y2_gar_res, y3_gar_res)
cop_res <- cop_res - min(cop_res)
cop_res <- cop_res / max(cop_res)
cop_fit1 <- RVineStructureSelect(as.copuladata(cop_res))

summary(cop_fit1)

# method by rescaling everything
copres <- cbind(y_gar_res, y2_gar_res, y3_gar_res, y4_gar_res, y5_gar_res, y6_gar_res)
copres <- copres - min(copres)
copres <- copres / max(copres)
cop_fit2 <- RVineStructureSelect(as.copuladata(copres))

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

# take one day of data
y11o <- pghyp(y[1:24], object=g11$best.model)
y22o <- pghyp(y2[1:24], object=g22$best.model)
y33o <- pghyp(y3[1:24], object=g33$best.model)
y44o <- pghyp(y4[1:24], object=g44$best.model)
y55o <- pghyp(y5[1:24], object=g55$best.model)
y66o <- pghyp(y6[1:24], object=g66$best.model)

svinedatao <- cbind(y11o, y22o, y33o, y44o, y55o, y66o)
colnames(svinedatao) <- c("M", "Tam", "Tal", "J", "O", "FM")
s_fit1 <- svinecop(svinedatao, p=7)
plot(s_fit1, var_names="use")

# 1 month of data
y11 <- pghyp(y[1:672], object=g11$best.model)
y22 <- pghyp(y2[1:672], object=g22$best.model)
y33 <- pghyp(y3[1:672], object=g33$best.model)
y44 <- pghyp(y4[1:672], object=g44$best.model)
y55 <- pghyp(y5[1:672], object=g55$best.model)
y66 <- pghyp(y6[1:672], object=g66$best.model)

svinedata <- cbind(y11, y22, y33, y44, y55, y66)
colnames(svinedata) <- c("M", "Tam", "Tal", "J", "O", "FM")
s_fit2 <- svinecop(svinedata, p=6)
plot(s_fit2, var_names="use")


# using monthly data
g11m <- stepAIC.ghyp(ym, silent=TRUE)
g22m <- stepAIC.ghyp(ym2, silent=TRUE)
g33m <- stepAIC.ghyp(ym3, silent=TRUE)
g44m <- stepAIC.ghyp(ym4, silent=TRUE)
g55m <- stepAIC.ghyp(ym5, silent=TRUE)
g66m <- stepAIC.ghyp(ym6, silent=TRUE)

y11m <- pghyp(ym[1:120], object=g11m$best.model)
y11m[is.na(y11m)] <- mean(y11m[!is.na(y11m)])
y22m <- pghyp(ym2[1:120], object=g22m$best.model)
y22m[is.na(y22m)] <- mean(y22m[!is.na(y22m)])
y33m <- pghyp(ym3[1:120], object=g33m$best.model)
y33m[is.na(y33m)] <- mean(y33m[!is.na(y33m)])
y44m <- pghyp(ym4[1:120], object=g44m$best.model)
y44m[is.na(y44m)] <- mean(y44m[!is.na(y44m)])
y55m <- pghyp(ym5[1:120], object=g55m$best.model)
y55m[is.na(y55m)] <- mean(y55m[!is.na(y55m)])
y66m <- pghyp(ym6[1:120], object=g66m$best.model)
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

y11mt <- pghyp(ymt[1:120], object=g11m$best.model)
y11mt[is.na(y11mt)] <- mean(y11mt[!is.na(y11mt)])
y22mt <- pghyp(ymt2[1:120], object=g22m$best.model)
y22mt[is.na(y22mt)] <- mean(y22mt[!is.na(y22mt)])
y33mt <- pghyp(ymt3[1:120], object=g33m$best.model)
y33mt[is.na(y33mt)] <- mean(y33mt[!is.na(y33mt)])
y44mt <- pghyp(ymt4[1:120], object=g44m$best.model)
y44mt[is.na(y44mt)] <- mean(y44mt[!is.na(y44mt)])
y55mt <- pghyp(ymt5[1:120], object=g55m$best.model)
y55mt[is.na(y55mt)] <- mean(y55mt[!is.na(y55mt)])
y66mt <- pghyp(ymt6[1:120], object=g66m$best.model)
y66mt[is.na(y66mt)] <- mean(y66mt[!is.na(y66mt)])

svinedatamt <- cbind(y11m, y22m, y33m, y44m, y55m, y66m, y11mt, y22mt, y33mt, y44mt, y55mt, y66mt)
colnames(svinedatamt) <- c("M W", "Tam W", "Tal W", "J W", "O W", "FM W", "M P", "Tam P", "Tal P", "J P", "O P", "FM P")
s_fitmt <- svinecop(svinedatamt, p=3)
plot(s_fitmt, var_names="use")



# ggplots to show the dependence 
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