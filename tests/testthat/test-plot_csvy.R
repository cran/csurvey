test_that("plot.csvy returns ggplot object for 1D example", {
  data("nhdat2", package = "csurvey")
  # specify the design
  dstrat <- survey::svydesign(ids = ~id, strata = ~str,
    data = nhdat2, weight = ~wt)
  
  set.seed(1)
  ans <- csvy(chol ~ incr(age), design = dstrat, n.mix = 10)
  
  p_default <- plot(ans)
  p_both <- plot(ans, type = "both")
  p_uncon <- plot(ans, type = "unconstrained")
  
  expect_s3_class(p_default, "ggplot")
  expect_s3_class(p_both,    "ggplot")
  expect_s3_class(p_uncon,   "ggplot")
})

test_that("plot.csvy handles two-factor layout and custom control", {
  # example with a binomial response
  library(NHANES)
  library(survey)
  data(NHANES)
  
  nh <- subset(NHANES, !is.na(Education) & !is.na(BMI) & !is.na(Weight))
  nh$DiabetesBin <- as.integer(nh$Diabetes == "Yes")
  nh$BMIgroup <- cut(nh$BMI, breaks = c(0, 18.5, 25, 30, 35, 40, Inf), labels = FALSE)
  
  # specify the design
  dsgn <- svydesign(ids = ~1, strata = ~BMIgroup, weights = ~Weight, data = nh)
  
  ans <- csvy(DiabetesBin ~ decr(Education) * incr(BMIgroup), design = dsgn,
              family = quasibinomial(link='logit'), n.mix=10)
  
  summary(ans)
  
  # first plot: BMIgroup on x1, Education on x2 
  p1 <- plot(ans, x1 = "BMIgroup", x2 = "Education")
  expect_s3_class(p1, "ggplot")
  
  # second plot with custom control
  ctl <- plot_csvy_control(x1size = 1.5, x2size = 2,
    angle  = 45, hjust  = 0.3)
  
  p2 <- plot(ans, x1 = "Education", x2 = "BMIgroup", control = ctl)
  expect_s3_class(p2, "ggplot")
})
