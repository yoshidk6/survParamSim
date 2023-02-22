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
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  # Lognormal distribution
  fit.lung.ln <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.ln <- surv_param_sim(fit.lung.ln, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.ln.raw <- extract_sim(sim.ln) %>% tibble::tibble()

  # Confirm that, at the median survival time predicted from survreg's predict function,
  # the created survival function would return 50% survival
  pred.median.surv <- predict(fit.lung.ln, newdata = newdata[1,])
  survfun.ln.1subj <- create_survfun(sim.ln.raw$lp[1], sim.ln$scale.vec[1])
  expect_equal(survfun.ln.1subj(pred.median.surv), 0.5)
  ## Make sure it also works with a vector of lp
  survfun.ln.1subj.rep <- create_survfun(sim.ln.raw$lp, sim.ln$scale.vec[1])
  expect_equal(survfun.ln.1subj.rep(pred.median.surv), 0.5)

})

#
# aaa <- extract_sim(sim.ln) %>%
#   tibble::tibble() %>%
#   dplyr::mutate(sex = factor(sex+1)) %>%
#   dplyr::mutate(trt = as.numeric(sex))
#
# sim.ln.nested <-
#   extract_sim(sim.ln) %>%
#   dplyr::mutate(sex = factor(sex)) %>%
#   dplyr::group_by(rep) %>%
#   tidyr::nest()

test_that("check grouping works", {

})



