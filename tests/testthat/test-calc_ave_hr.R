context("test-calc_ave_hr")

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


## >2 levels in treatment
newdata.3trt <-
  newdata %>%
  dplyr::mutate(trt = rep(c("B", "A", "CC"), length.out = nrow(.)),
                trtfct = factor(trt, levels = c("B", "A", "CC")))
# This will also test calculating grouping based on the variables not included in the model formula
fit.lung.3trt <- survreg(Surv(time, status) ~ sex + ph.ecog, data = newdata.3trt, dist = "lognormal")

sim.3trt <- surv_param_sim(fit.lung.3trt, newdata.3trt, n.rep, censor.dur)



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



test_that("calc_ave_hr_pi behavior check with manually calling", {
  # Manual calculation
  sim.raw <- extract_sim(sim) %>% tibble::tibble()
  data.each <- sim.raw %>%
    dplyr::filter(rep == 2) %>%
    dplyr::mutate(sex = forcats::fct_rev(factor(sex)))
  lp.vec.control <- data.each$lp[as.numeric(data.each$sex) == 1]
  lp.vec.treatment <- data.each$lp[as.numeric(data.each$sex) == 2]

  hr.manual.calc <-
    calc_ave_hr_from_lp(lp.vec.control, lp.vec.treatment, scale = sim$scale.vec[[2]],
                        dist = "lognormal", time.max = 1000)

  # Use `calc_ave_hr_pi()` function
  hr.pi <- calc_ave_hr_pi(sim, trt = "sex")

  # Expect the HR to be the same
  expect_equal(hr.manual.calc, 1/extract_hr(hr.pi)$HR[[2]])

  # Why it's not working..
  hr.pi <- calc_ave_hr_pi(sim, trt = "sex", trt.assign = "reverse")
})


test_that(">=3 trt groups", {
  sim.raw <- extract_sim(sim.3trt) %>% tibble::tibble()
  sim.nested <-
    sim.raw %>%
    dplyr::group_by(rep, trtfct) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::left_join(sim.3trt$scale.bs.df, by = "rep")

  df.lp.extracted <-
    sim.nested %>%
    dplyr::mutate(lp = purrr::map(data, function(x) x$lp)) %>%
    dplyr::select(-data)

  df.lp.control <-
    df.lp.extracted %>%
    dplyr::filter(as.numeric(trtfct) == 1) %>%
    dplyr::select(-trtfct)
  df.lp.treatment <-
    df.lp.extracted %>%
    dplyr::filter(as.numeric(trtfct) != 1) #%>%
    # dplyr::left_join(df.lp.control, by = c("rep")) # Need grouping variable

  dist <- sim.3trt$survreg$dist
  time.max <- 1000
  ## Scale

  df.surv.pdf.fun.control <-
    df.lp.control %>%
    dplyr::mutate(survfun.control = purrr::map2(lp, scale, function(x, y) create_survfun(lpvec = x, scale = y, dist = dist)),
                  pdf.control     = purrr::map2(lp, scale, function(x, y) create_pdf(lpvec = x, scale = y, dist = dist))) %>%
    dplyr::select(-lp, -scale)
  df.surv.pdf.fun.treatment <-
    df.lp.treatment %>%
    dplyr::mutate(survfun.trt = purrr::map2(lp, scale, function(x, y) create_survfun(lpvec = x, scale = y, dist = dist)),
                  pdf.trt     = purrr::map2(lp, scale, function(x, y) create_pdf(lpvec = x, scale = y, dist = dist))) %>%
    dplyr::select(-lp, -scale)

  ############## Need grouping variable ##########
  df.surv.pdf.fun.join <-
    dplyr::left_join(df.surv.pdf.fun.treatment, df.surv.pdf.fun.control, by = c("rep"))

  df.ahr <-
    df.surv.pdf.fun.join %>%
    dplyr::mutate(integrand1 = purrr::map2(survfun.control, pdf.trt, function(x, y) function(t){x(t) * y(t)}),
                  integrand2 = purrr::map2(survfun.trt, pdf.control, function(x, y) function(t){x(t) * y(t)})) %>%
    dplyr::mutate(term1 = purrr::map_dbl(integrand1, function(x) integrate(x, lower = 0, upper = time.max)$value),
                  term2 = purrr::map_dbl(integrand2, function(x) integrate(x, lower = 0, upper = time.max)$value)) %>%
    dplyr::mutate(HR = term1 / term2) %>%
    dplyr::select(rep, trtfct, HR)



  # Manual calculation
  sim.raw.rep2 <-
    sim.raw %>%
    dplyr::filter(rep == 2)

  lp.vec.control <- sim.raw.rep2$lp[as.numeric(sim.raw.rep2[["trtfct"]]) == 1]
  lp.vec.treatment <- sim.raw.rep2$lp[as.numeric(sim.raw.rep2[["trtfct"]]) == 2]

  survfun.control   <- create_survfun(lp.vec.control,   sim.3trt$scale.bs.df$scale[[2]], dist = dist)
  survfun.treatment <- create_survfun(lp.vec.treatment, sim.3trt$scale.bs.df$scale[[2]], dist = dist)
  pdf.control   <- create_pdf(lp.vec.control,   sim.3trt$scale.bs.df$scale[[2]], dist = dist)
  pdf.treatment <- create_pdf(lp.vec.treatment, sim.3trt$scale.bs.df$scale[[2]], dist = dist)

  integrand1 <- function(x){survfun.control(x) * pdf.treatment(x)}
  integrand2 <- function(x){survfun.treatment(x) * pdf.control(x)}

  term1 <- integrate(integrand1, lower=0, upper=time.max)$value
  term2 <- integrate(integrand2, lower=0, upper=time.max)$value

  ahr <- term1 / term2#
#   # Add scale
#
#   # First level in the factor serves as the reference group
#   lp.vec.control <- x$lp[as.numeric(x[[trt]]) == 1]
#   lp.vec.treatment <- x$lp[as.numeric(x[[trt]]) == 2]
#
#   ahr <-
#     calc_ave_hr_from_lp(lp.vec.control, lp.vec.treatment, scale = y,
#                         dist = "lognormal", time.max = time.max)

})

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



test_that("check trt.assign == reverse is behaving as expected", {

})


