context("test-calc_ave_hr")

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


test_that("check create_survfun function works as intended", {

  # Lognormal distribution
  fit.lung.ln <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.ln <- surv_param_sim(fit.lung.ln, newdata, n.rep = 1, censor.dur, coef.var = FALSE)

  # Confirm that, when using median survival time predicted from survreg's predict function,
  # the created survival function would return 50% survival
  pred.median.surv <- predict(fit.lung.ln, newdata = head(newdata, 1))
  survfun.ln.1subj <- create_survfun(sim.ln$lp.matrix[1,1], sim.ln$scale.vec[1])
  expect_equal(survfun.ln.1subj(pred.median.surv), 0.5)
  ## Make sure it also works with a vector of lp
  survfun.ln.1subj <- create_survfun(rep(sim.ln$lp.matrix[1,1], 5), sim.ln$scale.vec[1])
  expect_equal(survfun.ln.1subj(pred.median.surv), 0.5)

})
