
context("SDE class")

test_that("SDE constructor works", {
    n <- 100
    # expect_error(..., NA) checks that there is *no* error
    expect_error(SDE$new(formulas = list(mu = ~x, sigma = ~x), 
                         data = data.frame(ID = "ts1",
                                           x = rnorm(n), 
                                           z = rnorm(n), 
                                           time = 1:n), 
                         type = "BM", 
                         response = "z"), 
                 NA)
})

test_that("SDE fails if wrong data column names", {
    n <- 100
    f <- list(mu = ~x, sigma = ~x)
    
    # 1. Warning if ID is not in data
    data1 <- data.frame(x = rnorm(100),
                        z = rnorm(100),
                        time = 1:100)
    expect_warning(SDE$new(formulas = f, data = data1, 
                           type = "BM", response = "z"))
    
    # 2. Error if response is not in data
    data2 <- data.frame(ID = "ts1",
                        x = rnorm(100),
                        y = rnorm(100),
                        time = 1:100)
    expect_error(SDE$new(formulas = f, data = data2, 
                         type = "BM", response = "z"))
    
    # 3. Error if covariate is not in data
    data3 <- data.frame(ID = "ts1",
                        y = rnorm(100),
                        z = rnorm(100),
                        time = 1:100)
    expect_error(SDE$new(formulas = f, data = data3, 
                         type = "BM", response = "z"))
    
    # 4. Error if time is not in data
    data4 <- data.frame(ID = "ts1",
                        x = rnorm(100),
                        z = rnorm(100),
                        t = 1:100)
    expect_error(SDE$new(formulas = f, data = data4, 
                         type = "BM", response = "z"))
})

test_that("SDE members are correct size", {
    n <- 100
    data <- data.frame(ID = rep(paste0("ts", 1:(n/10)), each = 10),
                       x1 = rnorm(n),
                       x2 = rnorm(n),
                       z = rnorm(n),
                       time = 1:n)
    formulas <- list(mu = ~s(x1, k = 5, bs = "ts") + x2,
                     sigma = ~s(ID, bs = "re") + s(x2, k = 5, bs = "ts"))
    
    sde <- SDE$new(formulas = formulas, 
                   data = data, 
                   type = "BM",
                   response = "z")
    
    expect_equal(length(sde$coeff_fe()), 3)
    expect_equal(length(sde$coeff_re()), 18)
    expect_equal(length(sde$lambda()), 3)
    expect_equal(length(sde$vcomp()), 3)
})
