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

n.resample = c(10, 30)

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
               c(n.rep, 2))
})

