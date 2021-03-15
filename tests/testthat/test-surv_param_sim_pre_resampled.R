context("test-surv_param_sim_pre_resampled")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3) %>%
  tidyr::drop_na()
censor.dur <- c(200, 1100)

n.resample = c(100, 150)

sim.noresample <- surv_param_sim(object, newdata, n.rep, censor.dur)
sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = "sex")

newdata.resampled <-
  sim.resample$newdata.nona.sim %>%
  dplyr::select(-subj.sim, -n.resample)


sim.pre.resampled <-
  surv_param_sim_pre_resampled(object, newdata.orig = newdata, newdata.resampled = newdata.resampled, censor.dur = censor.dur)

test_that("Check N consistent across reps", {
  expect_error(surv_param_sim_pre_resampled(object,
                                            newdata.orig = newdata,
                                            newdata.resampled = dplyr::slice(newdata.resampled, -1L),
                                            censor.dur = censor.dur))
})

test_that("Error if NA present in covariates", {
  expect_error(surv_param_sim_pre_resampled(object,
                                            newdata.orig = tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)),
                                            newdata.resampled = newdata.resampled,
                                            censor.dur = censor.dur))
})

hr.pi <- calc_hr_pi(sim.pre.resampled, trt = "sex", group = "ph.ecog")
km.pi <- calc_km_pi(sim.pre.resampled, trt = "sex", group = "ph.ecog")

km.pi <- calc_km_pi(sim.pre.resampled, trt = "sex")
plot_km_pi(km.pi)

test_that("Error if generating hr_pi or km_pi with group variables not stratified for resampling", {
})


test_that("Check pre-resampled sim matches with surv_param_sim_resample results", {
})


# test_that("Extracted HR size matches", {
#   hr.pi <- calc_hr_pi(sim.resample, trt = "sex")
#
#   expect_equal(dim(extract_hr(hr.pi)),
#                c(n.rep, 2))
# })

