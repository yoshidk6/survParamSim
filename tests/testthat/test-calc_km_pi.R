context("test-calc_km_pi")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))
censor.dur <- c(200, 1100)


sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))
km.pi <- calc_km_pi(sim, group = "sex")

test_that("setting last obs time for simulation", {
  expect_equal(km.pi$t.last, 1022)
})

test_that("observed median time per group", {
  expect_equal(km.pi$obs.median.time$median, c(270, 426))
})

test_that("Extract median quantile in wide format", {
  median.pi.quantile <- extract_median_surv(km.pi, outtype ="wide")
  expect_equal(dim(median.pi.quantile), c(2, 6))
  expect_equal(names(median.pi.quantile), c("pi_low", "pi_med", "pi_high", "obs", "sex", "n"))

})

test_that("summary of prediction intervals for median survival time", {
  expect_equal(summary(km.pi)$median,
               c(213, 281, 331, 270, 298, 414, 535, 426),
               tolerance = 1)
  expect_equal(summary(km.pi)$n, rep(c(137, 90), each = 4))
})


test_that("predicted KM and median time per group", {
  km.pi$sim.median.time %>%
    dplyr::filter(rep == 1) %>%
    dplyr::pull(median) %>%
    expect_equal(c(295, 439), tolerance = 1)

  km.quantile <-
    km.pi$sim.km.quantile %>%
    dplyr::group_by(sex) %>%
    dplyr::slice(10) %>%
    dplyr::ungroup() %>%
    dplyr::select(pi_low, pi_high) %>%
    as.data.frame()

  km.quantile %>%
    dplyr::mutate(pi_low  = unname(pi_low),
                  pi_high = unname(pi_high)) %>%
    expect_equal(data.frame(pi_low  = c(0.777, 0.841),
                            pi_high = c(0.922, 0.951)),
                 tolerance = 0.01)
})



test_that("grouping and trt", {
  plot.km.sex <- plot_km_pi(calc_km_pi(sim, trt = "sex"))
  expect_doppelganger("km plot with sex as trt", plot.km.sex)

  plot.km.sex.ecog <- plot_km_pi(calc_km_pi(sim, group = c("sex", "ph.ecog")))
  expect_doppelganger("km plot with sex and ph.ecog as group", plot.km.sex.ecog)
})


test_that("long simulation time", {
  km.pi.longsim <- calc_km_pi(sim, group = "sex", simtimelast = 2000)

  expect_doppelganger("long sim time, not truncating with censor", plot_km_pi(km.pi.longsim, trunc.sim.censor = FALSE))
  expect_doppelganger("long sim time, truncating with censor", plot_km_pi(km.pi.longsim))
})



test_that("no group or trt", {
  km.pi.no.group.trt <- calc_km_pi(sim)

  summary(km.pi.no.group.trt) %>%
    dplyr::pull(median) %>%
    expect_equal(c(279, 320, 363, 310), tolerance = 1)

})


test_that("not calculating observed KM", {
  expect_doppelganger("no observed KM curves", plot_km_pi(calc_km_pi(sim, group = "sex", calc.obs = FALSE)))

})


