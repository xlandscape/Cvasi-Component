## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "../doc/figures/howto-",
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(cvasi)
library(dplyr)
#options(conflicts.policy=FALSE,warn.conflicts=FALSE)

## -----------------------------------------------------------------------------
# Create a new scenario object
myscenario <- Lemna_Schmitt()

# Get model name
get_model(myscenario)

# Set a custom tag to identify the scenario
myscenario %>%
  set_tag("Lab experiment #1") -> myscenario

# Get custom tag
get_tag(myscenario)

# The tag is also displayed when printing scenario properties
myscenario

# Accessing scenario slots and their default values
myscenario@forcings.req # forcings required for effect calculations
myscenario@endpoints    # available effect endpoints
myscenario@control.req  # are control runs required for effect calculation?

## ----eval=FALSE---------------------------------------------------------------
#  # Call the help page of effect scenarios class
#  ?scenarios

## ----eval=FALSE---------------------------------------------------------------
#  # Call the help page of the biomass transfer class
#  ?Transferable

## -----------------------------------------------------------------------------
# The example scenario `metsulfuron` based on the Lemna model by Schmitt et al. (2013)
# is modified by setting a new exposure series and initial state. Then, it is
# simulated.
  metsulfuron %>%
    set_noexposure() %>%  # set no exposure (i.e., a control run)
    set_init(c(BM = 50)) %>%  # set initial biomass
    simulate()

## -----------------------------------------------------------------------------
# Initialize the random numbers generator
set.seed(123)
# Generating a random exposure series spanning 14 days
random_conc <- runif(15, min=0, max=0.1)
exposure_profile <- data.frame(time=0:14, conc=random_conc)
# Run EPx calculations
minnow_it %>%
  set_exposure(exposure_profile) %>%  # set a specific exposure scenario
  epx()  # run EPx calculations

## -----------------------------------------------------------------------------
# Derive effect levels of all exposure windows
metsulfuron %>% 
  set_window(length=7, interval=1) %>%
  effect(max_only=FALSE)

## -----------------------------------------------------------------------------
# Only report effect levels larger than 1e-5 = 0.001%
metsulfuron %>% 
  set_window(length=7, interval=1) %>%
  effect(max_only=FALSE, marginal_effect=1e-5)

## -----------------------------------------------------------------------------
set.seed(123)
# Generate a random exposure series spanning 14 days
ts <- data.frame(time = 0:14,
                 conc = runif(15, min=0, max=0.1))

# Run EPx calculations for a window length of 7 days and a step size of 1 day
metsulfuron %>%
  set_exposure(ts) %>%
  set_window(length=7, interval=1) %>%
  effect(max_only=FALSE) -> results
results

# Create a plot of effects over time
library(ggplot2)
ggplot(results) +
  geom_point(aes(dat.start, BM*100)) +
  labs(x="Start of window (day)", y="Effect on biomass (%)")

## -----------------------------------------------------------------------------
metsulfuron %>%
  set_init(c(BM=1)) %>%
  set_noexposure() %>%
  set_transfer(interval=3, biomass=1) %>%
  simulate() -> result
result

library(ggplot2)
ggplot(result) +
  geom_line(aes(time, BM)) +
  labs(x="Time (days)", y="Biomass (g_dw/m2)", title="Biomass transfer every three days")


## -----------------------------------------------------------------------------
metsulfuron %>%
  set_init(c(BM=1)) %>%
  set_noexposure() %>%
  set_transfer(times=c(3, 6), biomass=c(1, 0.5)) %>%
  simulate() -> result2

library(ggplot2)
ggplot(result2) +
  geom_line(aes(time, BM)) +
  labs(x="Time (days)", y="Biomass (g_dw/m2)", title="Biomass transfer at custom time points")

## ----eval=FALSE---------------------------------------------------------------
#  # Call the help page of set_transfer
#  ?set_transfer

## ----warning=FALSE------------------------------------------------------------
library(dplyr)
# No exposure in control scenario
exposure <- data.frame(time=0:14, conc=0)

# set k_phot_fix, k_resp and Emax to the defaults recommended in Klein et al. (2022)
# for Tier 2C studies. Set initial biomass to 12.0 (original data is in fronds
# rather than biomass, but for sake of simplicity, no conversion was performed).
control <- metsulfuron %>% 
  set_param(c(k_phot_fix=TRUE, k_resp=0, Emax=1)) %>% 
  set_init(c(BM=12)) %>%
  set_exposure(exposure)
# `metsulfuron` is an example scenario based on Schmitt et al. (2013) and therefore
# uses the historical, superseded nomenclature by Schmitt et al. We recommend using
# the newer SETAC Lemna model by Klein et al. (2022) for current applications,
# see `Lemna_SETAC()` for details.

# Get control data from Schmitt et al. (2013)
obs_control <- Schmitt2013 %>%
  filter(ID == "T0") %>%
  select(t, BM=obs)

# Fit parameter `k_phot_max` to observed data: `k_phot_max` is a physiological
# parameter which is typically fitted to the control data
fit1 <- calibrate(
  x = control,
  par = c(k_phot_max = 1),
  data = obs_control,
  output = "BM",
  method="Brent", # Brent is recommended for one-dimensional optimization
  lower=0,        # lower parameter boundary
  upper=0.5       # upper parameter boundary
)
fit1$par

# Update the scenario with fitted parameter and simulate it
fitted_growth <- control %>% 
  set_param(fit1$par)
sim_mean <- fitted_growth %>%
  simulate() %>%
  mutate(trial="control")

# Exposure and observations in long format for plotting
treatments <- exposure %>%
  mutate(trial="control")
obs_mean <- obs_control %>%
  mutate(trial="control") %>%
  select(time=t, data=BM, trial)

# Plot results
plot_sd(
  model_base = fitted_growth,
  treatments = treatments,
  rs_mean = sim_mean,
  obs_mean = obs_mean
)

## -----------------------------------------------------------------------------
# get all the trials and exposure series from Schmitt et al. (2013) and create
# a list of calibration sets
Schmitt2013 %>%
  group_by(ID) %>%
  group_map(function(data, key) {
    exp <- data %>% select(t, conc)
    obs <- data %>% select(t, BM=obs)
    sc <- fitted_growth %>% set_exposure(exp) 
    caliset(sc, obs)
  }) -> cs

# Fit parameters on results of all trials at once
fit2 <- calibrate(
  cs,
  par=c(EC50=0.3, b=4.16, P_up=0.005),
  output="BM",
  method="L-BFGS-B",
  lower=c(0, 0.1, 0),
  upper=c(1000, 10, 0.1),
  verbose=FALSE
)
fit2$par

# Update the scenario with fitted parameter and simulate all trials
fitted_tktd <- fitted_growth %>%
  set_param(fit2$par)

treatments <- Schmitt2013 %>% select(time=t, conc, trial=ID)
rs_mean <- fitted_tktd %>%
  batch(treatments) %>%
  simulate()

# Observations in long format for plotting
obs_mean <- Schmitt2013 %>%
  select(time=t, data=obs, trial=ID)

# Plot results
plot_sd(
  model_base = fitted_tktd,
  treatments = treatments,
  rs_mean = rs_mean,
  obs_mean = obs_mean
)

## -----------------------------------------------------------------------------
# using the calibration set and calibrated parameters from the previous Lemna 
# Schmitt example, a likelihood profiling is done
# We set parameter boundaries to constrain the likelihood profiling within these 
# boundaries (this is optional)

# conduct profiling
res <- lik_profile(x = cs, # the calibration set
                   par = fit2$par, # the parameter values after calibration
                   output = "BM", # the observational output against which to 
                                  # compare the model fit
                   bounds = list(EC50_int = list(0.1, 4), 
                                 b = list(1, 5),
                                 P = list(0.0001, 0.2)))
# visualise profiling results
plot_lik_profile(res)

# access 95% confidence intervals of profiled parameters
res$EC50$confidence_interval
res$b$confidence_interval


## ----eval=FALSE---------------------------------------------------------------
#  # Call the help page for more information about the parameter space explorer
#  ?explore_space
#  

## -----------------------------------------------------------------------------

# conduct space exploration
res_space <- explore_space(
  x = cs,
  par = fit2$par,
  res = res, # output of the likelihood profiling function
  output = "BM")

# visualize the parameter space
plot_param_space(res_space)


## -----------------------------------------------------------------------------
# A base scenario is created for the whole period, but parameters are only
# valid until day 7
sc1 <- metsulfuron %>%
  set_times(0:14)

# A parameter change occurs at day 7: k_loss increases from 0 to 0.1 d-1
sc2 <- sc1 %>%
  set_param(c(k_loss=0.1))

# Combine both scenarios to a sequence with a break at *t = 7*:
# scenario `sc1` will be simulated for *t = [0, 7]*, and
# scenario `sc2` will be simulated for *t = [7, 14]*.
sq <- sequence(list(sc1, sc2), breaks=7)
simulate(sq)

## ----eval=FALSE---------------------------------------------------------------
#  # Call the help page of `sequence`
#  ?sequence

## -----------------------------------------------------------------------------
# Simulations with a maximum solver step length of hmax=0.01
metsulfuron %>%
  set_times(0:7) %>%
  simulate(hmax=0.1)

## ----eval=FALSE---------------------------------------------------------------
#  future::plan(future::multisession)

## ----eval=FALSE---------------------------------------------------------------
#  vignette("future-1-overview", package="future")

## ----warning=FALSE------------------------------------------------------------
## An exemplary implementation of the GUTS-RED-SD TKTD model ##
# Model ODEs following the deSolve specification
sd_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dDw <- kd * (Cw(t) - Dw)            # dDw/dt, scaled damage
    dH <- kk * max(0, Dw - z)           # dH/dt, cumulative hazard
    
    list(c(dDw, dH))                    # return derivatives
  })
}

## ----eval=FALSE---------------------------------------------------------------
#  vignette("deSolve", package="deSolve")

## ----warning=FALSE------------------------------------------------------------
## Properties of a sample scenario ##
init <- c(Dw=0, H=0)                         # initial state
times <- 0:5                                 # output time points [0,5]
param <- c(kd=22, hb=0.01, z=0.5, kk=0.08)   # model parameters
exp <- data.frame(time=0:5,                  # exposure time-series
                  conc=c(0,1,1,0.5,0.6,0.2))

# Create a linear interpolation of the exposure time-series
expf <- approxfun(x=exp$time, y=exp$conc, method="linear", rule=2)
# Extend parameter set by interpolated exposure series
paramx <- as.list(c(param, Cw=expf))

# Numerically solve the ODE
deSolve::ode(y=init, times=times, parms=paramx, func=sd_ode)

## -----------------------------------------------------------------------------
## Integrate a new model class into the package workflow ##
# Create a unique class that derives from 'EffectScenario'
setClass("MyGuts", contains="EffectScenario")

# Create an object of the new class and assign scenario properties
new("MyGuts", name="My custom model") %>%
  set_init(init) %>%
  set_times(times) %>%
  set_param(param) %>%
  set_exposure(exp, reset_times=FALSE) %>%
  set_endpoints("L") -> myscenario

myscenario

## -----------------------------------------------------------------------------
# the actual function calling deSolve can have a different signature
solver_myguts <- function(scenario, ...) {
  # get relevant data from scenario
  init <- scenario@init
  param <- scenario@param
  times <- scenario@times
  exp <- scenario@exposure@series
  if(nrow(exp) == 1) { # extend exposure series to have at least two rows
    row2 <- exp[1,]
    row2[[1]] <- row2[[1]]+1
    exp <- rbind(exp, row2)
  }
  
  # Create a linear interpolation of the exposure time-series
  expf <- approxfun(x=exp[,1], y=exp[,2], method="linear", rule=2)
  # Extend parameter set by interpolated exposure series
  paramx <- as.list(c(param, Cw=expf))
  
  # pass everything to ODE solver for numerical integration
  as.data.frame(deSolve::ode(y=init, times=times, parms=paramx, func=sd_ode, ...))
}

## Overload the solver() function ##
# The functions signature, i.e. the number and names of its arguments, must stay as is
setMethod("solver", "MyGuts", function(scenario, ...) solver_myguts(scenario, ...))

## -----------------------------------------------------------------------------
myscenario %>% simulate()

## -----------------------------------------------------------------------------
## Overload effect endpoint calculation ##
# fx() is called by effect()
setMethod("fx", "MyGuts", function(scenario, ...) fx_myguts(scenario, ...))

# @param scenario Scenario object to assess
# @param window Start & end time of the moving exposure window
# @param ... any additional parameters
fx_myguts <- function(scenario, window, ...) {
  # simulate the scenario (it is already clipped to the moving exposure window)
  out <- simulate(scenario, ...)
  # only use model state at the end of the simulation
  out <- tail(out, 1)
  # calculate survival according to EFSA Scientific Opinion on TKTD models
  # p. 33, doi:10.2903/j.efsa.2018.5377
  t <- unname(window[2] - window[1])
  survival <- exp(-out$H) * exp(-scenario@param$hb * t)
  return(c("L"=survival))
}

# Derive effect levels for our sample scenario
myscenario %>% effect()

## -----------------------------------------------------------------------------
## Get all model parameters
# Lemna model parameters taken from file 'mm2.r' included
# in supporting material of Schmitt et al. (2013)
param_study <- list(
  #     - Effect -
  Emax     = 0.784,   # [same as conc. data]  maximum effect concentration
  EC50     = 0.3,     # [same as conc. data]  Midpoint of effect curve
  b        = 4.16,    # [-]  Slope of effect curve
  #     - Toxicokinetics -
  P_up     = 0.0054,  # [cm/d]  Permeability for uptake
  AperBM   = 1000,    # [cm^2/g_dw]  Frond area/dry weight
  Kbm      = 0.75,    # []  Biomass(fw)-water partition coefficient
  P_Temp   = F,       # Switch for temperature dependence of cuticle permeability
  MolWeight = 381,    # [g/mol] Molmass of molecule (determines Q10_permeability)
  #     - Fate of biomass -
  k_phot_fix  = F,    # If True k_G_max is not changed by environmental factors
  k_phot_max  = 0.47, # [1/d]  Maximum photosynthesis rate
  k_resp   = 0.05,    # [1/d]  Respiration rate at ref. temperature
  k_loss   = 0.0,     # [1/d]  Some rate of loss (e.g. Flow rate)
  #     - Temperature dependence -
  Tmin     = 8.0  ,   # [°C]  Minimum growth temperature 
  Tmax     = 40.5 ,   # [°C]  Maximum growth temperature
  Topt     = 26.7 ,   # [°C]  Optimum growth temperature 
  t_ref    = 25,      # [°C]  Reference temperature for respiration rate
  Q10      = 2,       # [-]   Temperature dependent factor for respiration rate
  #     - Light dependence (Linear) -
  k_0      = 3 ,      # [1/d]  Intercept of linear part
  a_k      = 5E-5 ,   # [(1/d)/(kJ/m^2/d)]   Slope of linear part
  #     - Phosphorus dependence (Hill like dependence) -
  C_P      = 0.3,     # [mg/L]  Phosporus concentration in water
  CP50     = 0.0043,  # [mg/L]  P-conc. where growth rate is half
  a_P      = 1,       # []  Hill coefficient
  KiP      = 101,     # [mg/L]  P-inhibition constant for very high P-conc.
  #     - Nitrogen dependence (Hill like dependence) -
  C_N      = 0.6,     # [mg/L]  Nitrogen concentration in water
  CN50     = 0.034,   # [mg/L]  N-conc. where growth rate is half
  a_N      = 1,       # []  Hill coefficient
  KiN      = 604,     # [mg/L]  n-inhibition constant for very high P-conc.
  #     - Density dependence -
  BM50     = 176,     # [g_dw/m^2] Cut off BM
  #     - Others -
  mass_per_frond = 0.0001,  # [g_dw/frond]  Dry weight per frond
  BMw2BMd  = 16.7     # [g_fresh/g_dry]  Fresh- / dryweight 
)


## Define the forcing variables
forc_temp <- data.frame(t = 0, tmp = 12) # [°C]  Current temperature (may also be a table)
forc_rad  <- data.frame(t = 0, rad = 15000) # [kJ/m²/d]  Radiation  (may also be given as table)


## Define a simple exposure pattern
# t 0..6  concentration 1 ug/L
# t 7..14 concentration 0 ug/L
exposure <- data.frame(time = 0:14, 
                       conc = c(rep(1, 7), rep(0, 8))
                       )


## Set initial values 
# given in file 'mmc2.r' of Schmitt et al. (2013)
init <- c(
  BM       = 50,     # [g_dw/m^2]  Dry Biomass dry weight per m2
  E        = 1,      # (0-1)  (Toxic) Effect = Factor on growth rate  (Range: 0 - 1, 1=no effect)
  M_int    = 0       # [ug]   Amount of toxicant in biomass
)

## create a scenario object, containing the model (with parameters) and the exposure time-series
Lemna_Schmitt() %>%               # the Lemna model by Schmitt et al. (2013)
  set_tag("metsulfuron") %>%      # set a tage for the specific implementation of the model
  set_init(init) %>%              # set the starting values (as prepared above)
  set_param(param_study) %>%      # set the parameters (as prepared above)
  set_exposure(exposure) %>%      # set the exposure scenario (exposure time series)
  set_forcings(temp=forc_temp, rad=forc_rad) -> metsulfuron # set the external forcing 
# variables, and save everything as an object



## ----out.width = "100%"-------------------------------------------------------
## simulate with model, under a range of different exposure scenarios
# create several exposure scenarios
exp_scen <- data.frame(time = Schmitt2013$t,
                       conc = Schmitt2013$conc,
                       trial = Schmitt2013$ID)
# simulate for all these scenarios
results <- metsulfuron %>%
  batch(exp_scen) %>%
  simulate()

# plot results
plot_sd(
  model_base = metsulfuron,
  treatments = exp_scen,
  rs_mean = results,
)


## simulate with model, under a range of different exposure scenarios, and including 
## a biomass transfer
# simulate for all scenarios
results <- metsulfuron %>%
  set_transfer(interval = c(5), biomass = 10) %>% # implement a biomass transfer every 5 days
  batch(exp_scen) %>%
  simulate()

# plot results
plot_sd(
  model_base = metsulfuron,
  treatments = exp_scen,
  rs_mean = results,
)


