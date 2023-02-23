context("test-calc_ave_hr_pi")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")

object <- fit.lung
n.rep  <-  30

## ph.ecog == 3 only has one subject
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::filter(ph.ecog != 3)
censor.dur <- c(200, 1100)


sim <- surv_param_sim(object, newdata, n.rep, censor.dur)
hr.pi <- calc_ave_hr_pi(sim, trt = "sex", simtimelast = 1000, boot.subj = FALSE)


## >2 levels in treatment
newdata.3trt.2 <-
  newdata %>%
  dplyr::mutate(trt = rep(c("B", "A", "CC"), length.out = nrow(.)),
                trtfct = factor(trt, levels = c("B", "A", "CC")))
# This will also test calculating grouping based on the variables not included in the model formula
fit.lung.3trt.2 <- survreg(Surv(time, status) ~ sex + ph.ecog, data = newdata.3trt.2, dist = "lognormal")

sim.3trt.2 <- surv_param_sim(fit.lung.3trt.2, newdata.3trt.2, n.rep, censor.dur)



test_that("check create_survfun function works as intended", {
  newdata.1subj.rep <- dplyr::bind_rows(newdata[1,], newdata[1,], newdata[1,])

  # Lognormal distribution
  fit.lung.ln <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung, dist = "lognormal")
  sim.ln <- surv_param_sim(fit.lung.ln, newdata.1subj.rep, n.rep = 1, censor.dur, coef.var = FALSE)
  sim.ln.raw <- extract_sim(sim.ln) %>% tibble::tibble()

  # Confirm that, at the median survival time predicted from survreg's predict function,
  # the created survival function would return 50% survival
  pred.median.surv <- predict(fit.lung.ln, newdata = newdata[1,])
  survfun.ln.1subj <- create_survfun(sim.ln.raw$lp[1], sim.ln$scale.bs.df$scale[1])
  expect_equal(survfun.ln.1subj(pred.median.surv), 0.5)
  ## Make sure it also works with a vector of lp
  survfun.ln.1subj.rep <- create_survfun(sim.ln.raw$lp, sim.ln$scale.bs.df$scale[1])
  expect_equal(survfun.ln.1subj.rep(pred.median.surv), 0.5)

})



test_that("calc_ave_hr_pi behavior check with manually calling", {
  # Manual calculation
  sim.raw <- extract_sim(sim) %>% tibble::tibble()
  data.each <- sim.raw %>%
    dplyr::filter(rep == 2) %>%
    dplyr::mutate(sex = forcats::fct_rev(factor(sex)))
  lp.vec.control <- data.each$lp[as.numeric(data.each$sex) == 1]
  lp.vec.treatment <- data.each$lp[as.numeric(data.each$sex) == 2]

  hr.manual.calc <-
    calc_ave_hr_from_lp(lp.vec.control, lp.vec.treatment, scale = sim$scale.bs.df$scale[[2]],
                        dist = "lognormal", simtimelast = 1000)

  # Use `calc_ave_hr_pi()` function
  hr.pi <- calc_ave_hr_pi(sim, trt = "sex", simtimelast = 1000, boot.subj = FALSE)
  expect_equal(hr.manual.calc, 1/extract_hr(hr.pi)$HR[[2]])

  hr.pi <- calc_ave_hr_pi(sim, trt = "sex", simtimelast = 1000, trt.assign = "reverse", boot.subj = FALSE)
  expect_equal(hr.manual.calc, extract_hr(hr.pi)$HR[[2]])
})


test_that(">=3 trt arms", {

  hr.pi <- calc_ave_hr_pi(sim.3trt.2, trt = "trtfct", simtimelast = 1000, boot.subj = FALSE)

  # Manual calculation
  sim.raw <- extract_sim(sim.3trt.2) %>% tibble::tibble()

  sim.raw.rep2 <-
    sim.raw %>%
    dplyr::filter(rep == 2)

  lp.vec.control <- sim.raw.rep2$lp[as.numeric(sim.raw.rep2[["trtfct"]]) == 1]
  lp.vec.trt1 <- sim.raw.rep2$lp[as.numeric(sim.raw.rep2[["trtfct"]]) == 2]
  lp.vec.trt2 <- sim.raw.rep2$lp[as.numeric(sim.raw.rep2[["trtfct"]]) == 3]

  survfun.control <- create_survfun(lp.vec.control, sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)
  survfun.trt1    <- create_survfun(lp.vec.trt1,    sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)
  survfun.trt2    <- create_survfun(lp.vec.trt2,    sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)
  pdf.control <- create_pdf(lp.vec.control,   sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)
  pdf.trt1    <- create_pdf(lp.vec.trt1,      sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)
  pdf.trt2    <- create_pdf(lp.vec.trt2,      sim.3trt.2$scale.bs.df$scale[[2]], dist = sim.3trt.2$survreg$dist)

  integrand1 <- function(x){survfun.control(x) * pdf.trt1(x)}
  integrand2 <- function(x){survfun.trt1(x) * pdf.control(x)}
  term1 <- integrate(integrand1, lower=0, upper=1000)$value
  term2 <- integrate(integrand2, lower=0, upper=1000)$value
  ahr1 <- term1 / term2
  integrand1 <- function(x){survfun.control(x) * pdf.trt2(x)}
  integrand2 <- function(x){survfun.trt2(x) * pdf.control(x)}
  term1 <- integrate(integrand1, lower=0, upper=1000)$value
  term2 <- integrate(integrand2, lower=0, upper=1000)$value
  ahr2 <- term1 / term2


  # Expect the HR to be the same
  expect_equal(extract_hr(hr.pi) %>% dplyr::filter(rep == 2) %>% dplyr::pull(HR),
               c(ahr1, ahr2))

})


test_that("check grouping works", {

  # Manual calculation
  sim.raw <- extract_sim(sim) %>% tibble::tibble()
  data.each <- sim.raw %>%
    dplyr::filter(rep == 2,
                  ph.ecog == 0) %>%
    dplyr::mutate(sex = forcats::fct_rev(factor(sex)))
  lp.vec.control <- data.each$lp[as.numeric(data.each$sex) == 1]
  lp.vec.treatment <- data.each$lp[as.numeric(data.each$sex) == 2]

  hr.manual.calc <-
    calc_ave_hr_from_lp(lp.vec.control, lp.vec.treatment, scale = sim$scale.bs.df$scale[[2]],
                        dist = "lognormal", simtimelast = 1000)

  # Use `calc_ave_hr_pi()` function
  hr.pi <- calc_ave_hr_pi(sim, trt = "sex", group = "ph.ecog", trt.assign = "reverse", simtimelast = 1000, boot.subj = FALSE)
  expect_equal(hr.manual.calc,
               extract_hr(hr.pi) %>%
                 dplyr::filter(rep == 2, ph.ecog == 0) %>%
                 dplyr::pull(HR))


})



test_that("Check quantile calculation", {
  hr.pi.quantile <- extract_hr_pi(hr.pi)
  hr.pi.raw <- extract_hr(hr.pi)

  expect_equal(hr.pi.raw %>%
                 dplyr::reframe(quantile = quantile(HR, probs = c(0.025, 0.5, 0.975))) %>%
                 dplyr::pull(quantile) %>% as.numeric(),
               hr.pi.quantile %>% dplyr::pull(HR) %>% .[1:3])
})


