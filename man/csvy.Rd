\name{csvy}
\alias{csvy}
\alias{barplot.csvy}
%barplot.svystat {survey}
\alias{coef.csvy}
\alias{confint.csvy}
\alias{deff.csvy}
%svyby {survey}
%\alias{dotchart.csvy}
\alias{ftable.csvy}
%ftable.svystat {survey}
\alias{SE.csvy}
%svyby {survey}
\alias{summary.csvy}
\alias{svycontrast.csvy}
%svycontrast {survey}
\alias{vcov.csvy}
%surveysummary {survey}
\alias{predict.csvy}
\alias{plotpersp.csvy}
\alias{fitted.csvy}

\title{Design-based estimation of domain means with monotonicity constraints.}

\description{
The csvyby function performs design-based domain mean estimation with monotonicity and block-monotone shape constraints.
}

\usage{
csvy(formula, design, subset = NULL, family = stats::gaussian(),
          nD = NULL, level = 0.95, n.mix = 100L, test = TRUE,...)
\method{barplot}{csvy}(height, beside = TRUE,...)
\method{coef}{csvy}(object,...)
\method{confint}{csvy}(object, parm, level = 0.95, type = c("link", "response"), ...)
\method{deff}{csvy}(object,...)
%\method{dotchart}{csvy}(x,...,pch = 19)
\method{ftable}{csvy}(x,...)
\method{SE}{csvy}(object,...)
\method{summary}{csvy}(object,...)
\method{svycontrast}{csvy}(stat, contrasts,...)
\method{vcov}{csvy}(object,...)
\method{predict}{csvy}(object, newdata, type = c("link", "response"), se.fit
                 = FALSE, level = 0.95,...)
\method{fitted}{csvy}(object,...)
%\method{plotpersp}{csvyby}(object,...)
}


\arguments{
  \item{formula}{
    A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. A predictor can be a non-parametrically modelled variable with a monotonicity or block ordering restriction, or a combination of both. In terms of a non-parametrically modelled predictor, the user is supposed to indicate the relationship between the domain mean and a predictor \eqn{x} in the following way:

    Assume that \eqn{\mu} is the vector of domain means and \eqn{x} is a predictor:
    \itemize{
    \item{incr(x): }{\eqn{\mu} is increasing in \eqn{x}.}
    \item{decr(x): }{\eqn{\mu} is decreasing in \eqn{x}.}
    \item{block.Ord(x): }{\eqn{\mu} is has a block ordering in \eqn{x}.}
    \item{conc(x): }{\eqn{\mu} is concave in \eqn{x}.}
    \item{conv(x): }{\eqn{\mu} is convex in \eqn{x}.}
    \item{incr.conc(x): }{\eqn{\mu} is increasing and concave in \eqn{x}.}
    \item{incr.conv(x): }{\eqn{\mu} is increasing and convex in \eqn{x}.}
    \item{decr.conc(x): }{\eqn{\mu} is decreasing and concave in \eqn{x}.}
    \item{decr.conv(x): }{\eqn{\mu} is decreasing and convex in \eqn{x}.}
    }
}
  \item{design}{
A survey design, which must be specified by the svydesign routine in the survey package.
}

\item{subset}{
  Expression to select a subpopulation.
}

  \item{nD}{
  Total number of domains.
}

\item{family}{
  A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. There are four families: Gaussian, binomial, poisson, and Gamma.
}

%\item{start}{TBA}

% \item{amat}{
%A \eqn{k \times M} matrix imposing shape constraints in each dimension, where \eqn{M} is the
%total number of domains. A subroutine in the csurvey package will create a constraint matrix according to shape constraints specified in the formula. The default is amat = NULL.
%}

  \item{level}{
Confidence level of the approximate confidence surfaces. The default is 0.95.
}
  \item{n.mix}{
 The number of simulations used to get the approximate confidence intervals or surfaces. If \code{n.mix = 0}, no simulation will be done and the face of the final projection will be used to compute the covariance matrix of the constrained estimate. The default is n.mix = 100L.
}
  \item{test}{
A logical scalar. If test == TRUE, then the p-value for the test \eqn{H_0:\theta} is in \eqn{V} versus \eqn{H_1:\theta} is in \eqn{C} is returned. \eqn{C} is the constraint cone of the form \eqn{\{\beta: A\beta \ge 0\}}, and \eqn{V} is the null space of \eqn{A}. The default is test = TRUE.
}

 \item{...}{
 This term includes two other arguments: deff and multicore. deff = TRUE will request a design effect from svymean. multicore = TRUE will use multicore package to distribute subsets over multiple processors.
}

The \code{coef} function returns estimated systematic component of a csvy object.

The \code{confint} function returns the confidence interval of a csvy object. If type = "response", then the interval is for the mean; if type = "link", then the interval is for the systematic component.

\item{parm}{
An argument in the generic confint function in the stats package. For now, this argument is not in use.
}

The following arguments are used in the \code{predict} function.

  \item{object}{A csvy object.}

  \item{newdata}{A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}

  \item{type}{
  If the response is Gaussian, \code{type = "response"} and \code{type = "link"} give the predicted mean; if the response is binomial, poisson or Gamma, \code{type = "response"} gives the predicted mean, and \code{type = "link"} gives the predicted systematic component.
}

\item{se.fit}{
Logical switch indicating if confidence intervals are required.
}

The following arguments are used in the \code{barplot} function. See \code{\link[survey]{barplot.svystat}} for more details.

\item{height}{Analysis result.}

\item{beside}{Grouped, rather than stacked, bars.}

%\item{object}{TBA}

%\item{parm}{TBA}

%\item{df}{TBA}

The following arguments are used in the \code{ftable} function. See \code{\link[survey]{ftable.svystat}} for more details.

\item{x}{A csvy object.}

%\item{pch}{TBA}

The following arguments are used in the \code{svycontrast} function. See \code{\link[survey]{svycontrast}} for more details.

\item{stat}{A csvy object.}

\item{contrasts}{A vector or list of vectors of coefficients, or a call or list of calls.}

}



\value{
The output is a list of values used for estimation, inference and visualization. Main output include:

\item{survey.design}{The survey design used in the model.}
\item{etahat}{Estimated shape-constrained domain systematic component.}
\item{etahatu}{Estimated unconstrained domain systematic component.}
\item{muhat}{Estimated shape-constrained domain means.}
\item{muhatu}{Estimated unconstrained domain means.}
\item{lwr}{Approximate lower confidence band or surface for the shape-constrained domain mean estimate.}
\item{upp}{Approximate upper confidence band or surface for the shape-constrained domain mean estimate.}
\item{lwru}{Approximate lower confidence band or surface for the unconstrained domain mean estimate.}
\item{uppu}{Approximate upper confidence band or surface for the unconstrained domain mean estimate.}
\item{amat}{The \eqn{k \times M} constraint matrix imposing shape constraints in each dimension, where \eqn{M} is the total number of domains.}
\item{grid}{A \eqn{M \times p} grid, where \eqn{p} is the total number of predictors or dimensions.}
\item{nd}{A vector of sample sizes in all domains.}
\item{Ds}{A vector of the number of domains in each dimension.}
\item{acov}{Constrained mixture covariance estimate of domain means.}
\item{cov.un}{Unconstrained covariance estimate of domain means.}
\item{CIC}{The cone information criterion proposed in Meyer(2013a). It uses the "null expected degrees of freedom" as a measure of the complexity of the model. See Meyer(2013a) for further details of cic.}
\item{CIC.un}{The cone information criterion for the unconstrained estimator.}
\item{zeros_ps}{Index of empty domain(s).}
\item{nd}{Sample size of each domain.}
\item{pval}{p-value of the one-sided test.}
\item{family}{The family parameter defined in the formula.}
\item{df.residual}{The observed degree of freedom for the residuals of a csvy fit.}
\item{df.null}{The degree of freedom for the null model of a csvy fit.}
\item{domain}{Index of each domain in the data set contained in the survey.design object.}
\item{null.deviance}{The deviance for the null model of a csvy fit.}
\item{deviance}{The residual deviance of a csvy fit.}
\item{ans.unc_cp}{A data frame including the grid which is the combination of domains in each predictor, the domain mean estimate, and the constrained standard error.}
}



\details{
In a one dimensional situation, we assume that \eqn{\bar{y}_{U_t}} are non-decreasing over \eqn{T} domains. If this monotonicity is not used in estimation, the population domain means can be estimated by the Horvitz-Thompson estimator or the Hajek estimator. To use the monotonicity information, this csvy function starts from the
Hajek estimates \eqn{\bar{y}_{S_t} = (\sum_{k\in S_t}y_k/\pi_k)/N_t} and the isotonic estimator \eqn{(\hat{\theta}_1,\ldots,\hat{\theta}_T)^T} minimizes the weighted sum of squared deviations from the sample domain means over the set of ordered vectors; that is, \eqn{\bold{\hat{\theta}}} is the minimizer of \eqn{(\tilde{\bold{y}}_{S} - \bold{\theta})^T \bold{W}_s (\tilde{\bold{y}}_{S} - \bold{\theta})} subject to \eqn{\bold{A\theta} \geq \bold{0}}, where \eqn{\bold{W}_S} is the diagonal matrix with elements \eqn{\hat{N}_1/\hat{N},\ldots,\hat{N}_D/\hat{N}}, and \eqn{\hat{N} = \sum_{t=1}^T \hat{N}_t} and \eqn{\bold{A}} is a \eqn{m\times T} constraint matrix imposing the monotonicity constraint.

Domains can also be formed from multiple covariates. In that case, a grid will be used to represent the domains. For example, if there are two predictors \eqn{x_1} and \eqn{x_2}, and \eqn{x_1} has values on \eqn{D_1} domains: \eqn{1,\ldots,D_1}, \eqn{x_2} has values on \eqn{D_2} domains: \eqn{1,\ldots,D_2}, then the domains formed by \eqn{x_1} and \eqn{x_2} will be a \eqn{D_1\times D_2} by \eqn{2} grid.

To get \eqn{100(1-\alpha)\%} approximate confidence intervals or surfaces for the domain means, we apply the method in Meyer, M. C. (2018). \eqn{\hat{p}_J} is the estimated probability that the projection of \eqn{y_s} onto \eqn{\cal C} lands on \eqn{\cal F_J}, and the \eqn{\hat{p}_J} values are obtained by simulating many normal random vectors with estimated domain means and covariance matrix \eqn{I}, where \eqn{I} is a \eqn{M \times M} matrix, and recording the resulting sets \eqn{J}.

The user needs to provide a survey design, which is specified by the svydesign function in the survey package, and also a data frame containing the response, predictor(s), domain variable, sampling weights, etc. So far, only stratified sampling design with simple random sampling without replacement (STSI) is considered in the examples in this package.

Note that when there might be any empty domain, the user must specify the total number of domains in the \eqn{nD} argument.

For binomial and Poisson families use \code{family=quasibinomial()} and \code{family=quasipoisson()} to avoid a warning about non-integer numbers of successes.  The `quasi' versions of the family objects give
the same point estimates and standard errors and do not give the warning.
}


\seealso{
  \code{\link[cgam]{plotpersp}}, to create a 3D Plot for a csvy Object with at least two predictors.

  \code{\link[cgam]{incr}}, to specify an increasing order constraint in a csvy formula.

  \code{\link[cgam]{decr}}, to specify a decreasing order constraint in a csvy formula.

  \code{\link[cgam]{conc}}, to specify a concave order constraint in a csvy formula.

  \code{\link[cgam]{conv}}, to specify a concave order constraint in a csvy formula.

  \code{\link[cgam]{incr.conc}}, to specify an increasing-concave order constraint in a csvy formula.

  \code{\link[cgam]{decr.conv}}, to specify an decreasing-convex order constraint in a csvy formula.

  \code{\link[cgam]{decr.conc}}, to specify an decreasing-concave order constraint in a csvy formula.

  \code{\link[cgam]{incr.conv}}, to specify an increasing-convex order constraint in a csvy formula.

  \code{\link{block.Ord}}, to specify a blocking ordering order constraint in a csvy formula.

  \code{\link[survey]{svyby}}, to compute survey statistics on subsets of a survey defined by factors.

  \code{\link[survey]{svymean}}, to compute means for data from complex surveys.

 % \code{\link{csvy}}, to perform design-based domain mean estimation with monotonicity and block-monotone shape constraints when the response is gaussian, binomial, poisson, etc.

  \code{\link[survey]{svyglm}}, to fit a generalised linear model to data from a complex survey design, with inverse-probability weighting and design-based standard errors.

}


\references{
 Xu, X. and Meyer, M. C. (2021) One-sided testing of population domain means in surveys.

  Oliva, C., Meyer, M. C., and Opsomer, J.D. (2020) Estimation and inference of domain means subject to qualitative constraints. \emph{Survey Methodology}

  Meyer, M. C. (2018) A Framework for Estimation and Inference in Generalized Additive Models with Shape and Order Restrictions. \emph{Statistical Science \bold{33(4)}} 595--614.

  Wu, J., Opsomer, J.D., and Meyer, M. C. (2016) Survey estimation of domain means that respect natural orderings. \emph{Canadian Journal of Statistics \bold{44(4)}} 431--444.

  Meyer, M. C. (2013a) Semi-parametric additive constrained regression.
\emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715.

  Lumley, T. (2004) Analysis of complex survey samples. \emph{Journal of Statistical Software \bold{9(1)}} 1--19.

}

\author{
  Xiyue Liao
}

\examples{
%\dontrun{
# Example 1: monotonic in one dimension
data(api)
mcat <- apipop$meals
for(i in 1:10){mcat[trunc(apipop$meals/10) + 1 == i] <- i}
mcat[mcat == 100] <- 10
mcat <- factor(mcat)
M <- 10 # total number of domains

nsp<-c(200, 200, 200)  ## sample sizes per stratum
es<-sample(apipop$snum[apipop$stype=='E'&!is.na(apipop$avg.ed)&!is.na(apipop$api00)],nsp[1])
ms<-sample(apipop$snum[apipop$stype=='M'&!is.na(apipop$avg.ed)&!is.na(apipop$api00)],nsp[2])
hs<-sample(apipop$snum[apipop$stype=='H'&!is.na(apipop$avg.ed)&!is.na(apipop$api00)],nsp[3])
sid<-c(es, ms, hs)

pw <- 1:6194*0 + 4421 / nsp[1]
pw[apipop$stype == 'M'] <- 1018 / nsp[2]
pw[apipop$stype == 'H'] <- 755 / nsp[3]

fpc <- 1:6194*0 + 4421
fpc[apipop$stype == 'M'] <- 1018
fpc[apipop$stype == 'H'] <- 755

strsamp <- cbind(apipop, mcat, pw, fpc)[sid, ]
dstrat <- svydesign(ids = ~snum, strata = ~stype, fpc = ~fpc, data = strsamp, weight = ~pw)
rds <- as.svrepdesign(dstrat, type = "JKn")

ansc1 <- csvy(api00 ~ decr(mcat), design = rds, family = gaussian(), nD = M)
# summary(ansc1)


# Example 2: unconstrained in x1 and increasing in x2 and x3
%\dontrun{
D1 <- 5
D2 <- 5
D3 <- 6
Ds <- c(D1, D2, D3)
M <- cumprod(Ds)[3] # total number of domains

x1vec <- 1:D1
x2vec <- 1:D2
x3vec <- 1:D3
grid <- expand.grid(x1vec, x2vec, x3vec)
N <- M*100*4
Ns <- rep(N/M, M)

mu.f <- function(x) {
  mus <- x[1]^(0.25) + 4*exp(0.5 + 2*x[2]) / (1 + exp(0.5 + 2*x[2])) + sqrt(1/4 + x[3])
  mus <- as.numeric(mus$Var1)
  return (mus)
}

mus <- mu.f(grid)

H <- 4
nh <- c(180, 360, 360, 540)
n <- sum(nh)
Nh <- rep(N/H, H)

#generate population
y <- NULL
z <- NULL

for(i in 1:M){
  Ni <- Ns[i]
  mui <- mus[i]
  ei <- rnorm(Ni, 0, sd = 1)
  yi <- mui + ei
  y <- c(y, yi)
  zi <- i/M + rnorm(Ni, mean = 0, sd = 1)
  z <- c(z, zi)
}

x1 <- rep(grid[,1], times = Ns)
x2 <- rep(grid[,2], times = Ns)
x3 <- rep(grid[,3], times = Ns)
domain <- rep(1:M, times = Ns)

cts <- quantile(z, probs = seq(0, 1, length = 5))
strata <- 1:N*0
strata[z >= cts[1] & z < cts[2]] <- 1
strata[z >= cts[2] & z < cts[3]] <- 2
strata[z >= cts[3] & z < cts[4]] <- 3
strata[z >= cts[4] & z <= cts[5]] <- 4
freq <- rep(N/(length(cts) - 1), n)

w0 <- Nh/nh
w <- 1:N*0
w[strata == 1] <- w0[1]
w[strata == 2] <- w0[2]
w[strata == 3] <- w0[3]
w[strata == 4] <- w0[4]

pop <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, domain = domain, strata = strata, w = w)
ssid <- stratsample(pop$strata, c("1" = nh[1], "2" = nh[2], "3" = nh[3], "4" = nh[4]))
sample.stsi <- pop[ssid, ,drop = FALSE]
ds <- svydesign(id = ~1, strata = ~strata, fpc = ~freq, weights = ~w, data = sample.stsi)

#domain means are increasing w.r.t x1, x2 and block monotonic in x3
ord <- c(1, 1, 2, 2, 3, 3)
ans <- csvy(y ~ incr(x1)*incr(x2)*block.Ord(x3, order=ord), design = ds, family = gaussian(),
nD = M, test = FALSE, n.mix = 0)

#3D plot of estimated domain means: x1 and x2
plotpersp(ans)

#3D plot of estimated domain means: x3 and x2
plotpersp(ans, x3, x2)

#3D plot of estimated domain means: x3 and x2 for each domain of x1
plotpersp(ans, x3, x2, categ = "x1")

#3D plot of estimated domain means: x3 and x2 for each domain of x1
plotpersp(ans, x3, x2, categ = "x1", NCOL = 3)

%}

}

\keyword{constrained regression}
