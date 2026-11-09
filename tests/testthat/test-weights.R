test_that("exact conditional weights match independent quadrature", {
  d <- toy_competing()
  d$time[1] <- .35
  d$status[1] <- 2L
  cuts <- c(0, .5, 1)
  cfit <- cens_pamm(~ x, d, "time", "status", tau = 1, cut_points = cuts)
  prepared <- SHpamm:::shp_validate_data(d, "time", "status", 1L, 1)$data
  ped <- SHpamm:::shp_make_ped(prepared, cuts, 1L, "sh")
  got <- SHpamm:::shp_cond_weights(ped, cfit, 1L, derivative = TRUE)
  rows <- which(ped$id == 1L)
  X <- predict(cfit, newdata = ped[rows, , drop = FALSE], type = "lpmatrix")
  lambda <- exp(drop(X %*% coef(cfit)))
  T <- .35
  integrand <- function(t) {
    ifelse(t <= T, 1, exp(-lambda[1] * (pmin(t, .5) - T) -
                             lambda[2] * pmax(t - .5, 0)))
  }
  expected <- c(integrate(integrand, 0, .5, rel.tol = 1e-10)$value / .5,
                integrate(integrand, .5, 1, rel.tol = 1e-10)$value / .5)
  expect_equal(got$weight[rows], expected, tolerance = 1e-8)
  expect_true(all(got$weight >= 0 & got$weight <= 1))

  j <- min(2L, length(coef(cfit)))
  h <- 1e-5
  hi <- lo <- cfit
  hi$coefficients[j] <- hi$coefficients[j] + h
  lo$coefficients[j] <- lo$coefficients[j] - h
  fd <- (SHpamm:::shp_cond_weights(ped, hi, 1L)$weight -
           SHpamm:::shp_cond_weights(ped, lo, 1L)$weight) / (2 * h)
  expect_lt(max(abs(fd - got$dw[, j])), 1e-6)
})

test_that("exact marginal weights average a step censoring survival", {
  d <- toy_competing()
  fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 5)
  w <- fit$ped_data$weight
  expect_true(all(is.finite(w) & w >= 0 & w <= 1))
  expect_equal(w[fit$ped_data$event == 1L],
               rep(1, sum(fit$ped_data$event == 1L)))
})
