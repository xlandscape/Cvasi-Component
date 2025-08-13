## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "../doc/figures/manual-",
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(cvasi)
library(purrr)
library(dplyr)

## -----------------------------------------------------------------------------
library(cvasi)

## Sample base R workflow ##
# Create a scenario object of model GUTS-RED-IT
my_it <- GUTS_RED_IT()
# Set model parameters
my_it <- set_param(my_it, c(kd=1.2, hb=0, alpha=9.2, beta=4.3))
# Print scenario details
my_it

## -----------------------------------------------------------------------------
## Example 'tidyr' workflow ##
# the pipeline (%>%) symbol passes results to the next statement
GUTS_RED_IT() %>%
  set_param(c(kd=1.2, hb=0, alpha=9.2, beta=4.3))

## ----eval=FALSE,include=FALSE-------------------------------------------------
#  install.packages("cvasi", dependencies=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("remotes", dependencies=TRUE)
#  remotes::install_github("cvasi-tktd/cvasi", dependencies=TRUE)

## -----------------------------------------------------------------------------
library(cvasi)

# Create a new GUTS-RED-IT scenario and set its model parameters
GUTS_RED_IT() %>%
  set_param(c(kd=1.2, hb=0, alpha=9.2, beta=4.3))

## -----------------------------------------------------------------------------
# Example GUTS-RED-IT scenario derived from an acute fish toxicity study
# of the fathead minnow and Chlorpyrifos (Geiger et al. 1988)
minnow_it %>%
  simulate()

## ----eval=FALSE---------------------------------------------------------------
#  # Access the package help on GUTS-RED type models
#  ?"GUTS-RED-models"

## -----------------------------------------------------------------------------
# Define an exposure time-series
myexposure <- data.frame(time=c(0, 1, 1.01, 5),
                         conc=c(10, 10, 0, 0))

# Create and parameterize a scenario
GUTS_RED_IT() %>%
  set_param(c(kd=1.2, hb=0, alpha=9.2, beta=4.3)) %>%
  set_exposure(myexposure) -> myscenario

## -----------------------------------------------------------------------------
# Print information about the scenario
myscenario


## -----------------------------------------------------------------------------
# Update the exposure time-series but keep former output time points
myscenario %>%
  set_noexposure()

## -----------------------------------------------------------------------------
# Selected model parameters
myparam <- c(p_M=3211, v=0.023, k_J=0.63)
# Initial body length L
myinit <- c(L=0.02)
# Constant non-zero exposure
myexposure <- data.frame(time=0, conc=1.72)

DEB_abj() %>%
  set_param(myparam) %>%
  set_init(myinit) %>%
  set_exposure(myexposure, reset_times=FALSE) %>%
  set_times(0:10) %>%              # Output times 0,1,2,...,10
  set_mode_of_action(4) %>%        # Method of Action #4 to be activated
  set_window(length=3, interval=1) # Using moving exposure windows of length 3 days

## -----------------------------------------------------------------------------
# Example scenario of the Lemna TKTD model
metsulfuron %>%
  set_times(0:7) %>%
  simulate()

## -----------------------------------------------------------------------------
# Same simulation period, but with a smaller step length of 0.1
metsulfuron %>%
  set_times(seq(0, 7, 0.1)) %>%
  simulate() %>%
  tail()

## -----------------------------------------------------------------------------
# Original simulation period, but with a maximum solver step length of hmax=0.01
metsulfuron %>%
  set_times(0:7) %>%
  simulate(hmax=0.01)

## ----eval=FALSE---------------------------------------------------------------
#  ?cvasi::simulate

## -----------------------------------------------------------------------------
# GUTS-RED-IT scenario of the fathead minnow and chlorpyrifos
minnow_it %>% effect()

## ----include=FALSE------------------------------------------------------------
# make sure that value in text are still up to date
testthat::expect_equal(minnow_it %>% effect() %>% dplyr::pull(L), 6.297e-5, tolerance=0.001)

## -----------------------------------------------------------------------------
# Setting up a custom scenario with a total simulation period of 14 days and
# an exposure window length of 7 to assess a trapezoidal exposure pattern
americamysis %>%
  set_window(7) %>%
  set_exposure(data.frame(t=c(0,3,4,7,8), c=c(0,0,3,3,0))) %>%
  set_times(0:14) -> mydeb

# Derive maximum effect level of all exposure windows
mydeb %>% effect()

## ----include=FALSE------------------------------------------------------------
# make sure that value in text are still up to date
testthat::expect_equal(mydeb %>% effect() %>% dplyr::pull(L), 0.0521, tolerance=0.001)

## -----------------------------------------------------------------------------
# Restrict assessed endpoints to structural length (L)
mydeb %>% 
  set_endpoints("L") %>%
  epx()

## ----include=FALSE------------------------------------------------------------
# make sure that value in text are still up to date
testthat::expect_equal(mydeb %>%  set_endpoints("L") %>% epx(ep_only=TRUE) %>% unlist(use.names=FALSE), c(1.162598, 1.711914), ignore_attr=TRUE, tolerance=0.01)

## -----------------------------------------------------------------------------
# Examine how the EP20 value is derived
minnow_it %>% epx(level=20, verbose=TRUE)

## -----------------------------------------------------------------------------
## Display selected scenario properties
metsulfuron@times  # Output times
metsulfuron@init   # Initial state

## Simulate the sample scenario
metsulfuron %>% simulate()

## -----------------------------------------------------------------------------
metsulfuron %>%
  simulate(nout=0) %>%
  head(5)

## -----------------------------------------------------------------------------
metsulfuron %>%
  simulate(nout=8) %>%
  head(5)

## -----------------------------------------------------------------------------
metsulfuron %>%
  set_times(0:7) %>%  # restrict scenario to the period [0,7]
  effect()

## -----------------------------------------------------------------------------
metsulfuron %>%
  set_window(length=7, interval=1) %>%   # enable moving exposure windows
  effect(max_only=FALSE)                 # return effects of all windows

## -----------------------------------------------------------------------------
metsulfuron %>%
  set_window(length=7, interval=1) %>%         # enable moving exposure windows
  effect(max_only=FALSE, marginal_effect=0.01) # return effects of all windows

## -----------------------------------------------------------------------------
metsulfuron %>% epx()

## -----------------------------------------------------------------------------
metsulfuron %>%
  set_endpoints("BM") %>%
  epx(verbose=TRUE)

