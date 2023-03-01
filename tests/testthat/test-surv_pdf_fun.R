context("test-surv_pdf_fun")

library(survival)
set.seed(12345)


## ph.ecog == 3 only has one subject
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3)
censor.dur <- c(200, 1100)





test_that("check create_survfun function works as intended", {
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  # Lognormal distribution
  fit.lung.ln <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.ln <- surv_param_sim(fit.lung.ln, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.ln.raw <- extract_sim(sim.ln) %>% tibble::tibble()

  # Confirm that, at the median survival time predicted from survreg's predict function,
  # the created survival function would return 50% survival
  ###### Also add different quantiles like 25th or 75th for additional validation #####
  pred.median.surv <- predict(fit.lung.ln, newdata = newdata[1,])
  survfun.ln.1subj <- create_survfun(sim.ln.raw$lp[1], sim.ln$scale.bs.df$scale[1])
  expect_equal(survfun.ln.1subj(pred.median.surv), 0.5)
  ## Make sure it also works with a vector of lp
  survfun.ln.1subj.rep <- create_survfun(sim.ln.raw$lp, sim.ln$scale.bs.df$scale[1])
  expect_equal(survfun.ln.1subj.rep(pred.median.surv), 0.5)

})

test_that("check create_survfun function works as intended", {

  t.out <- seq(100, 600, by = 100)
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  # Lognormal distribution
  fit.lung.ln <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.ln <- surv_param_sim(fit.lung.ln, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.ln.raw <- extract_sim(sim.ln) %>% tibble::tibble()

  surv.vec.manual <- plnorm(q=t.out, meanlog=sim.ln.raw$lp, sdlog=exp(sim.ln$scale.bs.df$scale[1]), lower=FALSE)
  survfun.ln.1subj.rep <- create_survfun(sim.ln.raw$lp, sim.ln$scale.bs.df$scale[1])
  expect_equal(survfun.ln.1subj.rep(t.out), surv.vec.manual)

})


