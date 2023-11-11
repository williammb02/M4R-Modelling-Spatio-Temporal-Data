miami_tp <- import_climate_data("tp", miami_index[1], miami_index[2])
tampa_tp <- import_climate_data("tp", tampa_index[1], tampa_index[2])
tallahassee_tp <- import_climate_data("tp", tallahassee_index[1], tallahassee_index[2])

miami_u10 <- import_climate_data("u10", miami_index[1], miami_index[2])
tampa_u10 <- import_climate_data("u10", tampa_index[1], tampa_index[2])
tallahassee_u10 <- import_climate_data("u10", tallahassee_index[1], tallahassee_index[2])

miami_v10 <- import_climate_data("v10", miami_index[1], miami_index[2])
tampa_v10 <- import_climate_data("v10", tampa_index[1], tampa_index[2])
tallahassee_v10 <- import_climate_data("v10", tallahassee_index[1], tallahassee_index[2])


# extremogram analysis
uvm <- cbind(miami_u10, miami_v10)
utm <- cbind(miami_u10, miami_tp)
t_mt <- cbind(miami_tp, tampa_tp)
# nice extremal dependence
# extremogram2(uvm, 0.95, 0.95, 100, 1)
# extremogram2(uvm, 0.05, 0.05, 100, 2)
# extremogram2(uvm, 0.05, 0.95, 100, 3)
extremogram2(uvm, 0.95, 0.05, 100, 4)


extremogram2(utm, 0.95, 0.95, 100, 1)
extremogram2(utm, 0.05, 0.05, 100, 2)
extremogram2(utm, 0.05, 0.95, 100, 3)
extremogram2(utm, 0.95, 0.05, 100, 4)

extremogram2(t_mt, 0.95, 0.95, 100, 1)
extremogram2(t_mt, 0.05, 0.05, 100, 2)
extremogram2(t_mt, 0.05, 0.95, 100, 3)
extremogram2(t_mt, 0.95, 0.05, 100, 4)

extremogram1(miami_u10, 0.95, 100, 1)
permfn1(miami_u10, 0.95, 50, 2, 1, 100)
extremogram1(miami_u10, 0.05, 100, 2)
# stronger pattern in v
extremogram1(miami_v10, 0.95, 100, 1)
extremogram1(miami_v10, 0.05, 100, 2)

