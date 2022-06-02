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


sim <- surv_param_sim(object, newdata, n.rep, censor.dur)
hr.pi <- calc_hr_pi(sim, trt = "sex")


## >2 levels in treatment
newdata.3trt <-
  newdata %>%
  dplyr::mutate(trt = rep(c("B", "A", "CC"), length.out = nrow(.)),
                trt = factor(trt, levels = c("B", "A", "CC")),
                trto= factor(trt, levels = c("B", "A", "CC"), ordered = T))
# This will also test calculating grouping based on the variables not included in the model formula
fit.lung.3trt <- survreg(Surv(time, status) ~ sex + ph.ecog, data = newdata.3trt)

sim.3trt <- surv_param_sim(fit.lung.3trt, newdata.3trt, n.rep, censor.dur)






# Check trt variable values

test_that("error if you have NA in treatment", {
  newdata.with.trt.na <-
    newdata %>%
    dplyr::mutate(trt = sex)

  newdata.with.trt.na[1, "trt"] <-  NA

  sim.tmp <- surv_param_sim(object, newdata.with.trt.na, n.rep, censor.dur)
  expect_error(calc_hr_pi(sim.tmp, trt = "trt"),
               "`trt` cannot has NA values")

})



test_that("detect unused trt levels", {
  newdata.with.trt.three.factor <-
    newdata %>%
    dplyr::mutate(trt = sex)

  newdata.with.trt.three.factor$trt <- factor(newdata.with.trt.three.factor$trt, levels = c("1", "2", "3"))

  sim.tmp <- surv_param_sim(object, newdata.with.trt.three.factor, n.rep, censor.dur)
  expect_warning(calc_hr_pi(sim.tmp, trt = "trt"), "`trt` variable is factor and has unused levels")

})


test_that("not all groups have both treatment arms", {
  newdata.wo.both.trt <-
    dplyr::tibble(time = 0,
                  status = 1,
                  sex = c(1, 1, 2, 1, 1),
                  ph.ecog = c(1, 1, 1, 2, 2))

  sim.tmp <- surv_param_sim(object, newdata.wo.both.trt, n.rep, censor.dur)

  expect_warning(calc_hr_pi(sim.tmp, trt = "sex", group = "ph.ecog"),
               "HR was not calculable in at least one subgroup for the observed data")

})


test_that("check HR calculation", {
  hr.pi.raw <- extract_hr(hr.pi)

  expect_equal(dim(hr.pi.raw), c(30, 5))
  expect_equal(hr.pi.raw$HR[[1]], 0.708, tolerance = .001)
  expect_equal(hr.pi.raw$p.value.logrank[[1]], 0.0340, tolerance = .001)


  hr.pi.raw.group <- extract_hr(calc_hr_pi(sim, trt = "sex", group = "ph.ecog", trt.assign = "rev"))
  expect_equal(dim(hr.pi.raw.group), c(90, 6))
  expect_equal(hr.pi.raw.group$HR[[1]], 1.3, tolerance = .01)
})

test_that("Extract HR quantile in wide format", {
  hr.pi.quantile <- extract_hr_pi(hr.pi, outtype ="wide")
  expect_equal(dim(hr.pi.quantile), c(1, 5))
  expect_equal(names(hr.pi.quantile), c("sex", "pi_low", "pi_med", "pi_high", "obs"))
})

test_that("check TRT levels assignment", {
  expect_equal(calc_hr_pi(sim, trt = "sex", group = "ph.ecog", trt.assign = "rev")$trt.levels,
               c("2", "1"))
})


test_that("check summary", {
  hr.pi.summary <- summary(hr.pi)

  expect_equal(dim(hr.pi.summary), c(4, 4))
  expect_equal(hr.pi.summary$HR, c(0.46, 0.624, 0.836, 0.596), tolerance = .001)

})

test_that("test trt with >2 levels", {

  hr.pi.3trt <- calc_hr_pi(sim.3trt, trt = "trt", group = "sex", trt.assign = "reverse")

  hr.pi.quantile <- extract_hr_pi(hr.pi.3trt)
  expect_equal(dim(hr.pi.quantile), c(16, 5))
  expect_equal(hr.pi.quantile$HR[[5]], 0.758, tolerance = .001)
  expect_equal(hr.pi.quantile$trt[[5]], factor("A", levels = c("B", "A", "CC")))

  expect_error(calc_hr_pi(sim.3trt, trt = "trto"), "`trt` cannot be an ordered factor")

})


