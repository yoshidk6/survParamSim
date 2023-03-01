context("test-surv_pdf_fun")

library(survival)
set.seed(12345)


## ph.ecog == 3 only has one subject
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3)
censor.dur <- c(200, 1100)

newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])



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


test_that("lognormal: predicted survival from create_survfun matches with predict.survreg", {
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  fit.lung.test.surv <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.test.surv <- surv_param_sim(fit.lung.test.surv, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

  # Confirm that, at the 10th/median/80th survival time predicted from survreg's predict function,
  # the created survival function would return corresponding survival
  pred.surv.qtiles <- predict(fit.lung.test.surv, newdata = newdata[1,], type = "quantile", p = c(0.1, 0.5, 0.8))
  survfun.test.surv.1subj <- create_survfun(sim.test.surv.raw$lp[1], sim.test.surv$scale.bs.df$scale[1])
  expect_equal(survfun.test.surv.1subj(pred.surv.qtiles), 1 - c(0.1, 0.5, 0.8))

  ## Make sure it also works with a vector of lp
  survfun.test.surv.1subj.rep <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.bs.df$scale[1])
  expect_equal(survfun.test.surv.1subj.rep(pred.surv.qtiles), 1 - c(0.1, 0.5, 0.8))

})

test_that("validate survival function and pdf with flexsurv", {

  t.out <- seq(100, 600, by = 100)

  fit.lung.test.surv <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  fit.lung.test.flexsurv <- flexsurv::flexsurvreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")

  pred.flexsurv.surv <-
    predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "survival", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred)
  pred.flexsurv.hazard <-
    predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "hazard", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred)
  pred.flexsurv.pdf <- pred.flexsurv.surv * pred.flexsurv.hazard

  sim.test.surv <- surv_param_sim(fit.lung.test.surv, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

  survfun.test <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.bs.df$scale[1])
  pdf.test <- create_pdf(sim.test.surv.raw$lp, sim.test.surv$scale.bs.df$scale[1])

  expect_equal(survfun.test(t.out), pred.flexsurv.surv)
  expect_equal(pdf.test(t.out), pred.flexsurv.pdf)

})


test_that("Predicted survival from create_survfun matches with predict.survreg", {

  test_surv_fun <- function(newdata, fit) {
    t.out <- seq(100, 600, by = 100)

    fit.lung.test.flexsurv <- flexsurv::flexsurvreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = fit$dist)

    pred.flexsurv.surv <-
      predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "survival", times = t.out) %>%
      dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred)
    pred.flexsurv.hazard <-
      predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "hazard", times = t.out) %>%
      dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred)
    pred.flexsurv.pdf <- pred.flexsurv.surv * pred.flexsurv.hazard

    sim.test.surv <- surv_param_sim(fit, newdata, n.rep = 1, censor.dur, coef.var = FALSE)
    sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

    survfun.test <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.bs.df$scale[1], dist = sim.test.surv$survreg$dist)
    pdf.test <- create_pdf(sim.test.surv.raw$lp, sim.test.surv$scale.bs.df$scale[1], dist = sim.test.surv$survreg$dist)

    expect_equal(survfun.test(t.out), pred.flexsurv.surv)
    expect_equal(pdf.test(t.out), pred.flexsurv.pdf)
  }

  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal"))
  # test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "gaussian"))
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "weibull"))
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "exponential"))
  expect_error(test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "loglogistic")))

})


