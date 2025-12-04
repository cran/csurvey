test_that("csvy runs an example with one predictor", {
  data(nhdat2, package = 'csurvey')
  #specify the design:
  dstrat <- svydesign(ids = ~id, strata = ~str, data = nhdat2, weight = ~wt)
  
  #uncomment to use parallel computing:
  #options(csurvey.multicore=TRUE)
  #mixture-variance-covariance matrix is simulated
  set.seed(1)
  # monotonic in 1 predictor
  ans <- csvy(chol ~ incr(age), design = dstrat, n.mix=10)
  
  expect_s3_class(ans, "csvy")
  expect_true(is.numeric(ans$muhat))

  s <- summary(ans)
  expect_s3_class(s, "summary.csvy") 
})

test_that("csvy runs 2D monotone example with additional unconstrained factor", {
  data("nhdat2", package = "csurvey")
  
  #specify the design:
  dstrat <- survey::svydesign(ids = ~id, strata = ~str,
    data   = nhdat2, weight = ~wt)
  
  # mnotone in age and wcat, unconstrained in icat (income),
  set.seed(1)
  ans <- csvy(chol ~ incr(age) * incr(wcat) * icat, design = dstrat,
    test = FALSE, n.mix = 10)
  
  expect_s3_class(ans, "csvy")
  
  # domain grid and means should match in length
  expect_equal(length(ans$muhat), nrow(ans$grid))
  
  # Ds should record the number of domains per dimension,
  # and use the expected predictor names
  expect_true(all(c("age", "wcat", "icat") %in% names(ans$Ds)))
  expect_true(all(ans$Ds > 0))
  
  # CIC values should be numeric
  expect_true(is.numeric(ans$CIC))
  expect_true(is.numeric(ans$CIC.un))
})

