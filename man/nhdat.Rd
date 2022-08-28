\name{nhdat}
\alias{nhdat}
\title{NHANES Data}

\description{
The NHANES study provides health data for a sample of the U.S. population.
There are N=1680 observations with complete records for cholesterol level, age, height, waist size, and gender for adults ages 21-40.
}

\usage{data("nhdat")}

\format{
  A data frame with 1680 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{ID vector}
    \item{\code{chol}}{Cholesterol level: 1 (>= 200 (mmol/L) and 0 (< 200 (mmol/L)}
    \item{\code{wcat}}{A factor of categorized waist size}
    \item{\code{gender}}{Gender}
    \item{\code{age}}{A factor of categorized age}
    \item{\code{wt}}{A vector of stratified sampling weight}
    \item{\code{str}}{A vector of strata variable in the stratified sampling design}
  }
}


\examples{
library(csurvey)
data(nhdat)

#specify a stratified design
dstrat <- svydesign(ids = ~ id,  strata = ~ str, data = nhdat,  weight = ~ wt)

#constrained estimate: domain mean of cholesterol level is increasing in age and waist size
#M is the total number of domains
M <- 160
ans <- csvy(chol ~ incr(age) * incr(wcat) * gender, design = dstrat,
  nD = M, family = quasibinomial(link = "logit"), n.mix = 0, test = FALSE)

plotpersp(ans, categ = "gender", type = "response", NCOL = 2)

}


\keyword{datasets}
