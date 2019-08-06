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

test_that("predicted KM and median time per group", {
  km.pi$sim.median.time %>%
    dplyr::filter(rep == 1) %>%
    dplyr::pull(median) %>%
    expect_equal(c(295, 439), tolerance = 1)

  km.pi$sim.km.quantile %>%
    dplyr::group_by(sex) %>%
    dplyr::slice(10) %>%
    dplyr::ungroup() %>%
    dplyr::select(pi_low, pi_high) %>%
    as.data.frame() %>%
    expect_equal(data.frame(pi_low  = c(0.777, 0.841),
                            pi_high = c(0.922, 0.951)),
                 tolerance = 0.01)
})


#### Want to add tests for graphics!!
plot_km_pi(calc_km_pi(sim, trt = "sex"))
plot_km_pi(calc_km_pi(sim, group = c("sex", "ph.ecog")))


test_that("long simulation time", {
  km.pi <- calc_km_pi(sim, group = "sex", simtimelast = 2000)
  plot_km_pi(km.pi, cut.sim.censor = FALSE)
  plot_km_pi(km.pi)
  #### Add tests!

  expect_equal(2 * 2, 4)
})



test_that("no group or trt case", {
  km.pi <- calc_km_pi(sim)

  summary(km.pi) %>%
    dplyr::pull(median) %>%
    expect_equal(c(279, 320, 363, 310), tolerance = 1)

})


test_that("not calculating observed KM", {
  sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))

  km.pi <- calc_km_pi(sim, group = "sex", calc.obs = FALSE)
  #### Add tests!
  plot_km_pi(km.pi)

  expect_equal(2 * 2, 4)
})


