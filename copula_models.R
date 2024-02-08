library(VineCopula)
library(network)
library(kdecopula)

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
cop_fit3 <- RVineStructureSelect(copres2)
cop_fit3c <- RVineStructureSelect(copres2, type = "CVine")

cop_pdf <- RVinePDF(copres2, cop_fit3)

# plots
contour(cop_fit3)
plot(cop_fit3)

contour(cop_fit3c)
plot(cop_fit3c)

# hypothesis test of two copulas
c(cop_fit3$AIC, cop_fit3$BIC) < c(cop_fit3c$AIC, cop_fit3c$BIC)

RVineVuongTest(copres2, cop_fit3c, cop_fit3)

# simulate
set.seed(2024) 

tests <- RVineSim(5832, cop_fit3)
hist(qghyp(tests[,1], object=g1$best.model), breaks=100)


# visualise location of all the cities in Florida
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

ggplot() +
  ggtitle("Map of Florida and Cities in the D-Vine Copula Model") +
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
