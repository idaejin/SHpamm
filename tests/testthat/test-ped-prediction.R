test_that("PED assigns boundary events once and retains competing events", {
  d <- data.frame(time = c(.5, .25, .75, 1),
                  status = c(1L, 2L, 0L, 0L), x = 1:4)
  prepared <- SHpamm:::shp_validate_data(d, "time", "status", 1L, 1)$data
  sh <- SHpamm:::shp_make_ped(prepared, c(0, .5, 1), 1L, "sh")
  ce <- SHpamm:::shp_make_ped(prepared, c(0, .5, 1), 1L, "censor")
  expect_equal(sum(sh$event), 1L)
  expect_equal(sh$interval[sh$event == 1L], factor(1L, levels = 1:2))
  expect_equal(sum(sh$id == 2L), 2L)
  expect_equal(sum(ce$event), 1L)
  expect_equal(sum(ce$event[ce$id == 4L]), 0L)
})

test_that("CIFs are coherent and V2/V3 covariances are finite", {
  d <- toy_competing()
  fit <- sh_pamm(~ x, d, "time", "status", tau = 1,
                 cut_points = 6, weights_from = "conditional", cens_formula = ~ x)
  p <- predict_cif(fit, data.frame(x = 0), times = seq(0, 1, by = .2))
  expect_equal(p$cif[1], 0)
  expect_true(all(diff(p$cif) >= 0))
  expect_true(all(p$lower <= p$cif & p$cif <= p$upper))
  for (m in c("V2", "V3")) {
    V <- vcov(fit, method = m)
    expect_true(all(is.finite(V)))
    expect_equal(V, t(V), tolerance = 1e-8)
  }
  expect_length(censor_support(fit), nrow(d))
  diagnostic <- sh_diagnostic(fit, "x")
  expect_s3_class(diagnostic, "htest")
  expect_true(is.finite(diagnostic$p.value))
})
