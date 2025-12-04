## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, fig.height = 5
)
library(csurvey)
library(survey)

## -----------------------------------------------------------------------------
library(csurvey)
data(nhdat2, package = "csurvey")
dstrat <- svydesign(ids = ~id,  strata = ~str, data = nhdat2,  weight = ~wt)

## -----------------------------------------------------------------------------
ans <- csvy(chol ~ incr(age), design = dstrat, n.mix = 100)

## -----------------------------------------------------------------------------
cat("CIC (constrained):", ans$CIC, "\n")

## -----------------------------------------------------------------------------
cat("CIC (unconstrained):", ans$CIC.un, "\n")

## -----------------------------------------------------------------------------
plot(ans, type = "both")

## -----------------------------------------------------------------------------
set.seed(1)
ans <- csvy(chol ~ incr(age)*incr(wcat)*icat, design = dstrat)

## -----------------------------------------------------------------------------
domains <- data.frame(age = c(24, 35), wcat = c(2, 4), icat = c(2, 3))
pans <- predict(ans, newdata = domains, se.fit = TRUE)
cat("Predicted values, confidence intervals and standard errors for specified domains:\n")
print (pans)

## -----------------------------------------------------------------------------
ctl <- list(x1lab = "waist", x2lab = "income", subtitle.size = 8)
plot(ans, x1 = "wcat", x2 = "icat", control = ctl, domains = domains)

## -----------------------------------------------------------------------------
plot(ans, x1 = "wcat", x2 = "icat", control = ctl, type="unconstrained")

## -----------------------------------------------------------------------------
data(nhdat, package = "csurvey")
dstrat <- svydesign(ids = ~ id,  strata = ~ str, data = nhdat,  weight = ~ wt)
set.seed(1)
ans <- csvy(chol ~ incr(age) * incr(wcat) * gender, design = dstrat,
            family = binomial(link = "logit"), test = TRUE)

## -----------------------------------------------------------------------------
summary(ans)

## -----------------------------------------------------------------------------
ctl <- list(angle = 0, x1size = 2, x2size = 2, x1lab = "waist", x2_labels = c("male", "female"),
  subtitle.size=6)
plot(ans, x1 = "wcat", x2 = "gender", type="both", control = ctl)

## ----eval = FALSE-------------------------------------------------------------
# url1 <- "[https://github.com/xliaosdsu/csurvey-data/raw/main/nscg19.rda](https://github.com/xliaosdsu/csurvey-data/raw/main/nscg19.rda)"
# url2 <- "[https://github.com/xliaosdsu/csurvey-data/raw/main/nscg19_2.rda](https://github.com/xliaosdsu/csurvey-data/raw/main/nscg19_2.rda)"
# 
# dest1 <- tempfile(fileext = ".rda")
# dest2 <- tempfile(fileext = ".rda")
# 
# download.file(url1, destfile = dest1, mode = "wb")
# download.file(url2, destfile = dest2, mode = "wb")
# 
# load(dest1)  # creates object nscg19
# load(dest2)  # creates object nscg19_2

