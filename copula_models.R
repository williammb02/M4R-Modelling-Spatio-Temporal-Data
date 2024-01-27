library(VineCopula)
library(network)

# fit a copula model to residuals at 3 locations

cop_res <- cbind(y_gar_res, y2_gar_res, y3_gar_res)
cop_res <- cop_res - min(cop_res)
cop_res <- cop_res / max(cop_res)
cop_fit1 <- RVineStructureSelect(as.copuladata(cop_res))

summary(cop_fit1)


copres <- cbind(y_gar_res, y2_gar_res, y3_gar_res, y4_gar_res, y5_gar_res, y6_gar_res)
copres <- copres - min(copres)
copres <- copres / max(copres)
cop_fit2 <- RVineStructureSelect(as.copuladata(copres))