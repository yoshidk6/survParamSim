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
km.pi <- calc_ave_km_pi(sim, trt = "sex", group = "ph.ecog", boot.subj = FALSE)
# plot_km_pi(km.pi)


## >2 levels in treatment
newdata.3trt.2 <-
  newdata %>%
  dplyr::mutate(trt = rep(c("B", "A", "CC"), length.out = nrow(.)),
                trtfct = factor(trt, levels = c("B", "A", "CC")))
# This will also test calculating grouping based on the variables not included in the model formula
fit.lung.3trt.2 <- survreg(Surv(time, status) ~ sex + ph.ecog, data = newdata.3trt.2, dist = "lognormal")

sim.3trt.2 <- surv_param_sim(fit.lung.3trt.2, newdata.3trt.2, n.rep, censor.dur)

km.pi.3trt.2 <- calc_ave_km_pi(sim.3trt.2, trt = "trtfct", group = "ph.ecog", boot.subj = FALSE)
# plot_km_pi(km.pi.3trt.2)

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

test_that("calc_ave_km_pi behavior check with manually calling", {

  newdata.2subj <-
    newdata %>%
    dplyr::group_by(sex) %>%
    dplyr::slice_head() %>%
    dplyr::ungroup()
  sim.2subj <- surv_param_sim(object, dplyr::bind_rows(newdata.2subj, newdata.2subj), n.rep, censor.dur,
                              coef.var = FALSE)

  t.out <- seq(0, sim.2subj$t.last.orig.new, length.out = 100)

  sim.raw <- extract_sim(sim.2subj) %>% tibble::tibble()
  data.each <- sim.raw %>%
    dplyr::filter(rep == 2)

  surv.vec.manual <-
    data.each %>%
    dplyr::mutate(survfun =
                    purrr::map(lp, function(x)
                      function(t) plnorm(q=t, meanlog=x, sdlog=exp(sim.2subj$scale.bs.df$scale[[2]]), lower=FALSE))) %>%
    dplyr::mutate(km = purrr::map(survfun, function(x) x(t.out))) %>%
    head(2) %>%
    tidyr::unnest(km) %>%
    dplyr::pull(km)

  km.pi <- calc_ave_km_pi(sim.2subj, group = "ph.ecog", trt = "sex", boot.subj = FALSE, calc.obs = FALSE)

  expect_equal(surv.vec.manual, extract_km_pi(km.pi)$pi_med %>% as.numeric())

})



