
# This script isn't meant to present a useful analysis, but to provide
# basic code to create a model with multiple terms. This can be helpful
# to investigate components of model objects as implemented in smoothSDE.

library(smoothSDE)

n_ID <- 10
n_by_ID <- 100
n <- n_ID * n_by_ID

data <- data.frame(ID = rep(paste0("ts", 1:n_ID), each = n_by_ID),
                   Z = rnorm(n),
                   x1 = runif(n),
                   x2 = runif(n),
                   x3 = sample(paste0("cat", 1:3), size = n, replace = TRUE),
                   time = 1:n)

formulas <- list(mu = ~ x1 + s(x1, k = 5, bs = "cr") + s(x2, k = 5, bs = "cr"),
                 sigma = ~x1 + s(x2, by = x3, k = 5, bs = "cr"))

sde <- SDE$new(formulas = formulas, data = data, 
               type = "BM", response = "Z")

mats <- sde$make_mat()
