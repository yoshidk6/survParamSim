context("test-calc_km_pi")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))
## ph.ecog == 3 only has one subject, also remove ph.ecog==NA
newdata2 <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3)
censor.dur <- c(200, 1100)
## Test for trt with more than 2 categories
newdata3 <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(!is.na(ph.ecog)) %>%
  dplyr::mutate(ph.ecog = factor(ph.ecog, levels = c(3, 1, 2, 0)))
fit.lung.3 <- survreg(Surv(time, status) ~ sex, data = newdata3)

sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))
km.pi <- calc_km_pi(sim, group = "sex")

sim.newdata2 <- suppressWarnings(surv_param_sim(object, newdata2, n.rep, censor.dur))
sim.newdata3 <- suppressWarnings(surv_param_sim(fit.lung.3, newdata3, n.rep, censor.dur))

test_that("setting last obs time for simulation", {
  expect_equal(km.pi$t.last, 1022)
})

test_that("observed median time per group", {
  expect_equal(km.pi$obs.median.time$median, c(270, 426))
})

test_that("Extract median quantile in wide format", {
  median.pi.quantile <- extract_medsurv_pi(km.pi, outtype ="wide")
  expect_equal(dim(median.pi.quantile), c(2, 6))
  expect_equal(names(median.pi.quantile), c("pi_low", "pi_med", "pi_high", "obs", "sex", "n"))
  expect_equal(median.pi.quantile$pi_med, c(281.2, 414.8), tolerance = 0.01)
})

test_that("`extract_median_surv()` still works", {
  withr::local_options(list(lifecycle_verbosity = "quiet"))
  median.pi.quantile <- extract_median_surv(km.pi)
  expect_equal(median.pi.quantile$median[1:2], c(213.2723, 281.2060), tolerance = 0.01)
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

  km.pi.sex.ecog <- suppressWarnings(calc_km_pi(sim, group = c("sex", "ph.ecog")))
  plot.km.sex.ecog <- plot_km_pi(km.pi.sex.ecog)
  expect_doppelganger("km plot with sex and ph.ecog as group", plot.km.sex.ecog)
})


test_that("median survival delta", {
  km.pi.ph.ecog <- suppressWarnings(calc_km_pi(sim.newdata3, trt = "ph.ecog", group = "sex", trt.assign = "reverse"))

  medsurv.delta <- extract_medsurv_delta(km.pi.ph.ecog)

  expect_equal(dim(medsurv.delta), c(180, 4))
  expect_equal(medsurv.delta[[2, "median_delta"]], 70,
               tolerance = 1)
})

test_that("median survival delta prediction interval", {
  km.pi.newdata2 <- calc_km_pi(sim.newdata2, trt = "sex", group = "ph.ecog")

  extract_medsurv_delta_pi(km.pi.newdata2, outtype = "wide") %>%
    dplyr::select(pi_low, pi_high) %>%
    as.data.frame() %>%
    dplyr::mutate(pi_low  = unname(pi_low),
                  pi_high = unname(pi_high)) %>%
    expect_equal(data.frame(pi_low  = c(-32, 32, -54),
                            pi_high = c(647, 325, 235)),
                 tolerance = 0.01)

})



test_that("no group or trt", {
  km.pi.no.group.trt <- calc_km_pi(sim)

  summary(km.pi.no.group.trt) %>%
    dplyr::pull(median) %>%
    expect_equal(c(279, 320, 363, 310), tolerance = 0.01)

})


test_that("not calculating observed KM", {
  expect_doppelganger("no observed KM curves", plot_km_pi(calc_km_pi(sim, group = "sex", calc.obs = FALSE)))

})

test_that("warning with immature median time simulations", {
  sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur = c(50, 100)))
  expect_warning(calc_km_pi(sim, group = "sex"))

})


