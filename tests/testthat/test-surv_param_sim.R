context("test-surv_param_sim")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))
censor.dur <- c(200, 1100)


sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))

# sim <- surv_param_sim(object, newdata, n.rep = 1, censor.dur)


test_that("have NA in dataset", {
  expect_warning(surv_param_sim(object, newdata, n.rep, censor.dur),
                 "Not all subjects in `newdata`")
})


