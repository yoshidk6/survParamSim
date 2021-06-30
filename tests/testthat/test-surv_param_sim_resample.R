context("test-surv_param_sim_resample")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  tidyr::drop_na()
censor.dur <- c(200, 1100)

n.resample = c(100, 150)

sim.noresample <- surv_param_sim(object, newdata, n.rep, censor.dur)
sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = "sex")
sim.resample.nostrat <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample = 100)


test_that("Simulated data fame size matches", {
  expect_equal(dim(sim.resample$sim),
               c(n.rep * sum(n.resample), 4))
})

test_that("Extracted HR size matches", {
  hr.pi <- calc_hr_pi(sim.resample, trt = "sex")

  expect_equal(dim(extract_hr(hr.pi)),
               c(n.rep, 5))
})

test_that("Warning with n per subgroup not consistent", {
  expect_warning(calc_km_pi(sim.resample, trt = "sex", group = "ph.ecog"))
})

test_that("Extracted median surv size matches", {
  km.pi <- suppressWarnings(calc_km_pi(sim.resample, trt = "sex", group = "ph.ecog"))

  expect_equal(dim(extract_medsurv(km.pi)),
               c(196, 5))
  expect_equal(dim(extract_medsurv_pi(km.pi)),
               c(nrow(dplyr::distinct(dplyr::select(newdata, sex, ph.ecog)))*4, 8))
})

test_that("Expect warning for unbalanced subjects due to NA", {
  newdata.withna <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))

  sim.resample.withna <- suppressWarnings(surv_param_sim_resample(object, newdata.withna, n.rep, censor.dur, n.resample, strat.resample = c("sex")))
  expect_warning(calc_km_pi(sim.resample.withna))
})

test_that("Make sure km and hr calc works", {
  hr.pi <- calc_hr_pi(sim.resample, trt = "sex")
  km.pi <- calc_km_pi(sim.resample, trt = "sex")

  expect_equal(extract_hr_pi(hr.pi)$HR[[1]], 0.328, tolerance = .001)
  extract_km_pi(km.pi)
  plot_hr_pi(hr.pi)
  plot_km_pi(km.pi)
})

test_that("Make sure km and hr calc works with group", {
  hr.pi <- suppressWarnings(calc_hr_pi(sim.resample, trt = "sex", group = "ph.ecog"))
  km.pi <- suppressWarnings(calc_km_pi(sim.resample, trt = "sex", group = "ph.ecog"))

  expect_equal(extract_hr_pi(hr.pi)$HR[[5]], 0.361, tolerance = .001)
  extract_km_pi(km.pi)
  plot_km_pi(km.pi)
  plot_hr_pi(hr.pi)
})



