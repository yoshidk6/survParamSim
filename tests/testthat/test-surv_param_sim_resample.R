context("test-surv_param_sim_resample")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>% tidyr::drop_na()
censor.dur <- c(200, 1100)

n.resample = c(10, 30)

sim.noresample <- surv_param_sim(object, newdata, n.rep, censor.dur)
sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = "sex")


test_that("#subjects simulated after resampling", {
  expect_equal(nrow(sim.resample$sim),
               n.rep * sum(n.resample))
})



