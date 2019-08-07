context("test-calc_hr_pi")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30

## ph.ecog == 3 only has one subject
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3)
censor.dur <- c(200, 1100)


sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))
hr.pi <- calc_hr_pi(sim, trt = "sex")



# Check trt variable values

test_that("error if you have NA in treatment", {
  newdata.with.trt.na <-
    newdata %>%
    dplyr::mutate(trt = sex)

  newdata.with.trt.na[1, "trt"] <-  NA

  sim.tmp <- suppressWarnings(surv_param_sim(object, newdata.with.trt.na, n.rep, censor.dur))
  expect_error(calc_hr_pi(sim.tmp, trt = "trt"),
               "`trt` cannot has NA values")

})


test_that("error if trt has more than 2 values", {
  newdata.with.three.trt <-
    newdata %>%
    dplyr::mutate(trt = sex)

  newdata.with.three.trt[1, "trt"] <-  3

  sim.tmp <- suppressWarnings(surv_param_sim(object, newdata.with.three.trt, n.rep, censor.dur))

  expect_error(calc_hr_pi(sim.tmp, trt = "trt"),
               "`trt` should contain exactly two unique values")

})

test_that("check if trt is factor with three or more levels", {
  newdata.with.trt.three.factor <-
    newdata %>%
    dplyr::mutate(trt = sex)

  newdata.with.trt.three.factor$trt <- factor(newdata.with.trt.three.factor$trt, levels = c("1", "2", "3"))

  sim.tmp <- suppressWarnings(surv_param_sim(object, newdata.with.trt.three.factor, n.rep, censor.dur))
  expect_error(calc_hr_pi(sim.tmp, trt = "trt"),
               "`trt` should have only two factor levels")

})


test_that("HR plots", {
  vdiffr::expect_doppelganger("HR plot without group", plot_hr_pi(calc_hr_pi(sim, trt = "sex")))
  vdiffr::expect_doppelganger("HR plot by ph.ecog", plot_hr_pi(calc_hr_pi(sim, trt = "sex", group = "ph.ecog", trt.assign = "rev")))
})



