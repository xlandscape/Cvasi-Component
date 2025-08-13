library(hdf5r)
library(cvasi)
library(pbapply)

h5file <- commandArgs(TRUE)[1]
pboptions(type = "timer")

f <- H5File$new(h5file, "r+")
temp_data <- f[["/envir/tmp"]][]
rad_data <- f[["/envir/irr"]][]
Lemna_Schmitt(
  c(
    Emax = f[["/param/Emax"]][],
    AperBM = f[["/param/AperBM"]][],
    Kbm = f[["/param/Kbm"]][],
    P_Temp = f[["/param/P_Temp"]][],
    MolWeight = f[["/param/MolWeight"]][],
    k_phot_fix = f[["/param/k_phot_fix"]][],
    k_phot_max = f[["/param/k_phot_max"]][],
    k_resp = f[["/param/k_resp"]][],
    k_loss = f[["/param/k_loss"]][],
    Tmin = f[["/param/Tmin"]][],
    Tmax = f[["/param/Tmax"]][],
    Topt = f[["/param/Topt"]][],
    t_ref = f[["/param/t_ref"]][],
    Q10 = f[["/param/Q10"]][],
    k_0 = f[["/param/k_0"]][],
    a_k = f[["/param/a_k"]][],
    C_P = f[["/param/C_P"]][],
    CP50 = f[["/param/CP50"]][],
    a_P = f[["/param/a_P"]][],
    KiP = f[["/param/KiP"]][],
    C_N = f[["/param/C_N"]][],
    CN50 = f[["/param/CN50"]][],
    a_N = f[["/param/a_N"]][],
    KiN = f[["/param/KiN"]][],
    BM50 = f[["/param/BM50"]][],
    mass_per_frond = f[["/param/mass_per_frond"]][],
    BMw2BMd = f[["/param/BMw2BMd"]][],
    EC50 = f[["/param/EC50"]][],
    b = f[["/param/b"]][],
    P_up = f[["/param/P_up"]][]
  ),
  c(BM = f[["/init/BM"]][], E = f[["/init/E"]][], M_int = f[["/init/M_int"]][])
) |>
  set_forcings(
    temp = data.frame(t = seq_along(temp_data), tmp = temp_data),
    rad = data.frame(t = seq_along(rad_data), rad = rad_data)
  ) |>
  set_times(0:(f[["/times"]][] - 1)) ->
  scenario
days <- 1:f[["/envir/conc"]]$dims[2] / 24
effects <- pblapply(1:f[["/envir/conc"]]$dims[1], function(reach) {
  concentrations <- f[["/envir/conc"]][reach,]
  scenario |>
    set_exposure(data.frame(time = days, conc = concentrations), reset_times = FALSE) ->
    scenario
  scenario |>
    simulate() ->
    result
  f[["/out/BM"]][reach,] <- result$BM
  f[["/out/E"]][reach,] <- result$E
  f[["/out/M_int"]][reach,] <- result$M_int
  f[["/out/C_int"]][reach,] <- result$C_int
  f[["/out/FrondNo"]][reach,] <- result$FrondNo
  scenario |>
    epx()
})
f[["/effects/BM.EP10"]][] <- sapply(effects, \(x) x$BM.EP10)
f[["/effects/r.EP10"]][] <- sapply(effects, \(x) x$r.EP10)
f[["/effects/BM.EP50"]][] <- sapply(effects, \(x) x$BM.EP50)
f[["/effects/r.EP50"]][] <- sapply(effects, \(x) x$r.EP50)
f$close()
