context("test-extract_raw_sim")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>% tidyr::drop_na()
censor.dur <- c(200, 1100)

n.resample = c(10, 30)


test_that("make sure extraction works for both with or without resampling", {
  sim <- surv_param_sim(object, newdata, n.rep, censor.dur)
  sim.raw <- extract_sim(sim)

  sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = "sex")
  sim.raw.resample <- extract_sim(sim.resample)

  expect_equal(dim(sim.raw), c(6810, 6))
  expect_equal(dim(sim.raw.resample), c(1200, 6))
})
