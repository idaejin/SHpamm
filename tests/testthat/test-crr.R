test_that("linear SH-PAMM association agrees with Fine-Gray regression", {
  skip_if_not_installed("cmprsk")
  set.seed(44)
  n <- 1300
  x <- rnorm(n)
  t1 <- rexp(n, .27 * exp(.45 * x))
  t2 <- rexp(n, .18)
  cens <- rexp(n, .12)
  d <- data.frame(
    time = pmin(t1, t2, cens, 5),
    status = ifelse(t1 <= pmin(t2, cens, 5), 1L,
                    ifelse(t2 <= pmin(t1, cens, 5), 2L, 0L)),
    x = x
  )
  fit <- sh_pamm(~ x, d, "time", "status", tau = 5, cut_points = 35)
  fg <- cmprsk::crr(d$time, d$status, as.matrix(d["x"]),
                    failcode = 1, cencode = 0)
  expect_equal(unname(stats::coef(fit)["x"]), unname(fg$coef[1]),
               tolerance = .02)
})
