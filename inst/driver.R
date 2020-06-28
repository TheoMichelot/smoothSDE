## Driver script

library(smoothSDE)

formulas <- list(mu = ~ s(x1, k = 10, bs = "cs"), sigma = ~1)

n <- 10
data <- data.frame(ID = 1,
                   Z1 = cumsum(rnorm(n)),
                   Z2 = cumsum(rnorm(n)),
                   x1 = cumsum(rnorm(n)),
                   x2 = cumsum(rnorm(n)))

type <- list(Z1 = "BM")

my_sde <- SDE$new(formulas = formulas, data = data, type = type)

mats <- my_sde$make_mats()