## Driver script

library(smoothSDE)

###################
## Simulate data ##
###################
## Drift function
true_mu <- function(x) {
    mu <- rep(0, length(x))
    ind1 <- which(x < 0.5)
    mu[ind1] <- sin(2 * pi * x[ind1] / 0.5)
    mu[-ind1] <- 2 * (plogis(25 * (x[-ind1] - 0.5)) - 0.5)
    return(mu)
}

## Simulation function for dZ_t = mu(x1) dt + dW_t
sim_SDE <- function(times, x1) {
    ## Number of points to simulate
    n <- length(times)
    ## Time intervals
    dt <- diff(times)
    ## Compute drift from covariates
    all_mu <- true_mu(x1)
    
    ## Generate normal increments for process Z
    dZ <- rnorm(n = n-1, mean = all_mu[-n]*dt, sd = sqrt(dt))
    ## Compute process Z
    Z <- cumsum(c(0, dZ))
    
    return(Z)
}

## Times of simulation
n_all <- 1e5
T_max <- 1e3
times_all <- seq(0, T_max, length = n_all)

## Simulate covariate (random walk, normalised to be in (0, 1))
x1_raw <- cumsum(rnorm(n_all))
x1_all <- (x1_raw - min(x1_raw))/(max(x1_raw) - min(x1_raw))

## Simulate data
Z_all <- sim_SDE(times = times_all, x1 = x1_all)

## Thin observations *regularly* (for estimation)
n_obs <- 1e3
ind_obs <- seq(1, n_all, length = n_obs)
Z <- Z_all[ind_obs]
x1 <- x1_all[ind_obs]
times <- times_all[ind_obs]

#################
## Setup model ##
#################
formulas <- list(mu = ~ s(x1, k = 10, bs = "cs"), sigma = ~1)
data <- data.frame(ID = 1, Z = Z, x1 = x1, time = times)
type <- "BM"

my_sde <- SDE$new(formulas = formulas, data = data, type = type)

mats <- my_sde$make_mats()

my_sde$fit(silent = FALSE)


my_sde$tmb_obj()$fn(my_sde$tmb_obj()$par)

#' TODO:
#'  - get estimated parameters to check that estimation works
#'  - write predict method to predict on grid
#'  - write plot method
#'  - think about best way to implement link functions
