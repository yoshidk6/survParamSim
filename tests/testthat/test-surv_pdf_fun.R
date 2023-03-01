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

  surv.vec.manual <- plnorm(q=t.out, meanlog=sim.ln.raw$lp, sdlog=exp(sim.ln$scale.ln.bs.df$scale.ln[1]), lower=FALSE)
  survfun.ln.1subj.rep <- create_survfun(sim.ln.raw$lp, sim.ln$scale.ln.bs.df$scale.ln[1])
  expect_equal(survfun.ln.1subj.rep(t.out), surv.vec.manual)

})


test_that("gaussian: predicted survival from create_survfun matches with predict.survreg", {
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  fit.lung.test.surv <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "gaussian")
  sim.test.surv <- surv_param_sim(fit.lung.test.surv, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

  # Confirm that, at the 10th/median/80th survival time predicted from survreg's predict function,
  # the created survival function would return corresponding survival
  pred.surv.qtiles <- predict(fit.lung.test.surv, newdata = newdata[1,], type = "quantile", p = c(0.1, 0.5, 0.8))
  survfun.test.surv.1subj <- create_survfun(sim.test.surv.raw$lp[1], sim.test.surv$scale.ln.bs.df$scale[1],
                                            dist = fit.lung.test.surv$dist)
  expect_equal(survfun.test.surv.1subj(pred.surv.qtiles), 1 - c(0.1, 0.5, 0.8))

  ## Make sure it also works with a vector of lp
  survfun.test.surv.1subj.rep <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale[1],
                                                dist = fit.lung.test.surv$dist)
  expect_equal(survfun.test.surv.1subj.rep(pred.surv.qtiles), 1 - c(0.1, 0.5, 0.8))

})


test_that("Predicted survival from create_survfun matches with predict.flexsurvreg, lognormal distribution", {

  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])
  fit <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")

  t.out <- seq(100, 600, by = 100)

  fit.lung.test.flexsurv <- flexsurv::flexsurvreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = fit$dist)

  pred.flexsurv.surv <-
    predict(fit.lung.test.flexsurv, newdata = newdata.1subj.rep[1,], type = "survival", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_survival)
  pred.flexsurv.hazard <-
    predict(fit.lung.test.flexsurv, newdata = newdata.1subj.rep[1,], type = "hazard", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_hazard)
  pred.flexsurv.pdf <- pred.flexsurv.surv * pred.flexsurv.hazard

  sim.test.surv <- surv_param_sim(fit, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

  survfun.test <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)
  pdf.test <- create_pdf(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)

  expect_equal(survfun.test(t.out), pred.flexsurv.surv)
  expect_equal(pdf.test(t.out), pred.flexsurv.pdf)

})
test_that("Predicted survival from create_survfun matches with predict.flexsurvreg, loglogistic distribution", {

  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])
  fit <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "loglogistic")

  t.out <- seq(100, 600, by = 100)

  fit.lung.test.flexsurv <- flexsurv::flexsurvreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "llogis")

  pred.flexsurv.lp <-
    predict(fit.lung.test.flexsurv, newdata = newdata.1subj.rep[1,], type = "lp") %>%
    dplyr::pull(.pred_link)
  pred.flexsurv.surv <-
    predict(fit.lung.test.flexsurv, newdata = newdata.1subj.rep[1,], type = "survival", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_survival)
  pred.flexsurv.hazard <-
    predict(fit.lung.test.flexsurv, newdata = newdata.1subj.rep[1,], type = "hazard", times = t.out) %>%
    dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_hazard)
  pred.flexsurv.pdf <- pred.flexsurv.surv * pred.flexsurv.hazard

  sim.test.surv <- surv_param_sim(fit, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

  # Need custom tolerance. Linear predictor is slightly different between survreg and flexsurvreg
  survfun.test <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)
  pdf.test <- create_pdf(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)
  expect_equal(survfun.test(t.out), pred.flexsurv.surv, tolerance = 1e-3)
  expect_equal(pdf.test(t.out), pred.flexsurv.pdf, tolerance = 1e-3)

  # Using estimates from flexsurvreg to make sure function itself is doing what it's expected to do
  survfun.test <- create_survfun(log(pred.flexsurv.lp), -fit.lung.test.flexsurv$coefficients[["shape"]], dist = sim.test.surv$survreg$dist)
  pdf.test <- create_pdf(log(pred.flexsurv.lp), -fit.lung.test.flexsurv$coefficients[["shape"]], dist = sim.test.surv$survreg$dist)
  expect_equal(survfun.test(t.out), pred.flexsurv.surv)
  expect_equal(pdf.test(t.out), pred.flexsurv.pdf)

})


test_that("Predicted survival from create_survfun matches with predict.flexsurvreg", {

  test_surv_fun <- function(newdata, fit) {
    t.out <- seq(100, 600, by = 100)
    dist.flexsurv <- ifelse(fit$dist == "loglogistic", "llogis", fit$dist)

    fit.lung.test.flexsurv <- flexsurv::flexsurvreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = dist.flexsurv)

    pred.flexsurv.surv <-
      predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "survival", times = t.out) %>%
      dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_survival)
    pred.flexsurv.hazard <-
      predict(fit.lung.test.flexsurv, newdata = newdata[1,], type = "hazard", times = t.out) %>%
      dplyr::pull(.pred) %>% .[[1]] %>% dplyr::pull(.pred_hazard)
    pred.flexsurv.pdf <- pred.flexsurv.surv * pred.flexsurv.hazard

    sim.test.surv <- surv_param_sim(fit, newdata, n.rep = 1, censor.dur, coef.var = FALSE)
    sim.test.surv.raw <- extract_sim(sim.test.surv) %>% tibble::tibble()

    survfun.test <- create_survfun(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)
    pdf.test <- create_pdf(sim.test.surv.raw$lp, sim.test.surv$scale.ln.bs.df$scale.ln[1], dist = sim.test.surv$survreg$dist)

    expect_equal(survfun.test(t.out), pred.flexsurv.surv)
    expect_equal(pdf.test(t.out), pred.flexsurv.pdf)
  }

  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal"))
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "weibull"))
  test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "exponential"))

  # Gaussian not in flexsurv
  # test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "gaussian"))
  # Need custom tolerance. Linear predictor is slightly different between survreg and flexsurvreg
  # test_surv_fun(newdata.1subj.rep, survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "loglogistic"))

})


