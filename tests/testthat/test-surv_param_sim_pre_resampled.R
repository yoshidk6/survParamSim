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

n.resample = c(10, 30)

sim.noresample <- surv_param_sim(object, newdata, n.rep, censor.dur)
sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = c("sex"))

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

test_that("Expect warning for unbalanced subjects due to NA", {
  newdata.withna <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))

  sim.resample.withna <- suppressWarnings(surv_param_sim_resample(object, newdata.withna, n.rep=100, censor.dur, n.resample, strat.resample = c("sex")))
  newdata.resampled.withna <-
    sim.resample.withna$newdata.nona.sim %>%
    dplyr::select(-subj.sim, -n.resample)

  sim.pre.resample.withna <-
    expect_warning(surv_param_sim_pre_resampled(object,
                                                newdata.orig = newdata.withna,
                                                newdata.resampled = newdata.resampled.withna,
                                                censor.dur = censor.dur))

  expect_warning(calc_km_pi(sim.pre.resample.withna))
})




test_that("Expect warning for missing obs and call calc_km_pi or calc_hr_pi", {

  sim.pre.resampled.no.obs <-
    surv_param_sim_pre_resampled(object, newdata.resampled = newdata.resampled, censor.dur = censor.dur)

  calc_km_pi(sim.pre.resampled.no.obs, calc.obs = FALSE)
  calc_hr_pi(sim.pre.resampled.no.obs, calc.obs = FALSE, trt = "sex")

  expect_warning(calc_km_pi(sim.pre.resampled.no.obs))
  expect_warning(calc_hr_pi(sim.pre.resampled.no.obs, trt = "sex"))
})




test_that("Check pre-resampled sim matches with surv_param_sim or surv_param_sim_resample results", {

  skip('Too heavy to run routinely')

  n.rep  <-  1000
  newdata <-
    tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
    dplyr::filter(ph.ecog != 3) %>%
    tidyr::drop_na()

  n.resample = c(1000, 1200)


  # Just repeated samples to compare with surv_param_sim
  sim <- surv_param_sim(object, newdata, n.rep = n.rep, censor.dur = censor.dur)

  newdata.repeat <-
    dplyr::tibble(rep = 1:n.rep) %>%
    tidyr::expand_grid(newdata)
  sim.pre.repeated <- surv_param_sim_pre_resampled(object, newdata.orig = newdata, newdata.resampled = newdata.repeat, censor.dur = censor.dur)


  km.pi.3 <- calc_km_pi(sim.pre.repeated, trt = "sex")
  km.pi.4 <- calc_km_pi(sim, trt = "sex")
  plot_km_pi(km.pi.3)
  plot_km_pi(km.pi.4)
  extract_medsurv_pi(km.pi.3)
  extract_medsurv_pi(km.pi.4)


  # Resampled samples
  sim.resample <- surv_param_sim_resample(object, newdata, n.rep, censor.dur, n.resample, strat.resample = c("sex"))

  newdata.resampled <-
    sim.resample$newdata.nona.sim %>%
    dplyr::select(-subj.sim, -n.resample)
  sim.pre.resampled <- surv_param_sim_pre_resampled(object, newdata.orig = newdata, newdata.resampled = newdata.resampled, censor.dur = censor.dur)


  km.pi.2 <- calc_km_pi(sim.pre.resampled, trt = "sex")
  plot_km_pi(km.pi.2)
  km.pi.1 <- calc_km_pi(sim.resample, trt = "sex")
  plot_km_pi(km.pi.1)
  extract_medsurv_pi(km.pi.1)
  extract_medsurv_pi(km.pi.2)

})


test_that("Make sure PI get narrower with coef.var = FALSE", {

  skip('Too heavy to run routinely')

  n.rep  <-  300
  newdata <-
    tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
    dplyr::filter(ph.ecog != 3) %>%
    tidyr::drop_na()

  n.resample = c(100, 120)


  # Just repeated samples to compare with surv_param_sim
  sim <- surv_param_sim(object, newdata, n.rep = n.rep, censor.dur = censor.dur)

  newdata.repeat <-
    dplyr::tibble(rep = 1:n.rep) %>%
    tidyr::expand_grid(newdata)
  sim.pre.repeated <- surv_param_sim_pre_resampled(object, newdata.orig = newdata, newdata.resampled = newdata.repeat, censor.dur = censor.dur)
  sim.pre.repeated.nocoefvar <- surv_param_sim_pre_resampled(object, newdata.orig = newdata, newdata.resampled = newdata.repeat, censor.dur = censor.dur,
                                                             coef.var = FALSE)

  plot_km_pi(calc_km_pi(sim.pre.repeated))
  plot_km_pi(calc_km_pi(sim.pre.repeated.nocoefvar))
})

