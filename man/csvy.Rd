\name{csvy}
\alias{csvy}
\alias{summary.csvy}
\alias{vcov.csvy}
\alias{predict.csvy}
\alias{coef.csvy}
\alias{confint.csvy}
\alias{plotpersp.csvy}
\alias{barplot.csvy}


\title{
Estimation of Domain Means with Monotonicity or Convexity Constraints
}

\description{
The csvy function performs design-based domain mean estimation with monotonicity and block-monotone shape constraints.

For example, in a one dimensional situation, we assume that \eqn{\bar{y}_{U_t}} are non-decreasing over \eqn{T} domains. If this monotonicity is not used in estimation,
the population domain means can be estimated by the Horvitz-Thompson estimator or the Hajek estimator. To use the monotonicity information, this csvy function starts from the
Hajek estimates \eqn{\bar{y}_{S_t} = (\sum_{k\in S_t}y_k/\pi_k)/N_t} and the isotonic estimator \eqn{(\hat{\theta}_1,\ldots,\hat{\theta}_T)^T} minimizes the weighted sum of squared deviations from the sample domain means over the set of ordered vectors; that is, \eqn{\bold{\hat{\theta}}} is the minimizer of \eqn{(\tilde{\bold{y}}_{S} - \bold{\theta})^T \bold{W}_s (\tilde{\bold{y}}_{S} - \bold{\theta})} subject to \eqn{\bold{A\theta} \geq \bold{0}}, where \eqn{\bold{W}_S} is the diagonal matrix with elements \eqn{\hat{N}_1/\hat{N},\ldots,\hat{N}_D/\hat{N}}, and \eqn{\hat{N} = \sum_{t=1}^T \hat{N}_t} and \eqn{\bold{A}} is a \eqn{m\times T} constraint matrix imposing the monotonicity constraint. 

Domains can also be formed from multiple covariates. In that case, a grid will be used to represent the domains. For example, if there are two predictors \eqn{x_1} and \eqn{x_2}, and \eqn{x_1} has values on \eqn{D_1} domains: \eqn{1,\ldots,D_1}, \eqn{x_2} has values on \eqn{D_2} domains: \eqn{1,\ldots,D_2}, then the domains formed by \eqn{x_1} and \eqn{x_2} will be a \eqn{D_1\times D_2} by \eqn{2} grid.

To get \eqn{100(1-\alpha)\%} approximate confidence intervals or surfaces for the domain means, we apply the method in Meyer, M. C. (2018). \eqn{\hat{p}_J} is the estimated probability that the projection of \eqn{y_s} onto \eqn{\cal C} lands on \eqn{\cal F_J}, and the \eqn{\hat{p}_J} values are obtained by simulating many normal random vectors with estimated domain means and covariance matrix \eqn{I}, where \eqn{I} is a \eqn{M \times M} matrix, and recording the resulting sets \eqn{J}.

The user needs to provide a survey design, which is specified by the svydesign function in the survey package, and also a data frame containing the response, predictor(s), domain variable, sampling weights, etc. So far, only stratified sampling design with simple random sampling without replacement (STSI) is considered in the examples in this package.

Note that when there is any empty domain, the user must specify the total number of domains in the \eqn{nD} argument. 

}

\usage{
csvy(formula, design, family=stats::gaussian(), multicore=getOption("csurvey.multicore"), 
  level=0.95, n.mix=100L, test=FALSE, subset=NULL)
\method{summary}{csvy}(object,...)
\method{vcov}{csvy}(object,...)
\method{coef}{csvy}(object,...)
%\method{fitted}{csvy}(object,...)
\method{confint}{csvy}(object, parm=NULL, level = 0.95, type = c("link", "response"),...)
%\method{plotpersp}{csvy}(object,x1 = NULL, x2 = NULL, x1nm = NULL, x2nm = NULL, 
%    data = NULL, ci = "none", transpose = FALSE, main = NULL, 
%    categ = NULL, categnm = NULL, surface = c("C", "U"), type = c("link", 
%        "response"), col = "white", cex.main = 0.8, xlab = NULL, 
%    ylab = NULL, zlab = NULL, zlim = NULL, box = TRUE, axes = TRUE, 
%    th = NULL, ltheta = NULL, ticktype = "detailed", nticks = 5, 
%    palette = NULL, NCOL = NULL,...)
\method{predict}{csvy}(object, newdata = NULL, type = c("link", "response"), 
    se.fit = TRUE, level = 0.95, n.mix = 100,...)
}

\arguments{
  \item{formula}{
    A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. A predictor can be a non-parametrically modelled variable with a monotonicity or convexity restriction, or a combination of both. In terms of a non-parametrically modelled predictor, the user is supposed to indicate the relationship between the domain mean and a predictor \eqn{x} in the following way:

    Assume that \eqn{\mu} is the vector of domain means and \eqn{x} is a predictor:
    \describe{
    \item{incr(x): }{\eqn{\mu} is increasing in \eqn{x}.}
    \item{decr(x): }{\eqn{\mu} is decreasing in \eqn{x}.}
    \item{block.Ord(x): }{\eqn{\mu} is has a block ordering in \eqn{x}.}
    }
}
\item{design}{
A survey design, which must be specified by the svydesign routine in the survey package.
}

\item{family}{
  A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. There are four families used in csvy: Gaussian, binomial, poisson, and Gamma. 
}

\item{multicore}{
A parameter retrieving the current global option for "csurvey.multicore" and assigns its value to the multicore variable, allowing the function to respect a user-defined setting for parallel processing behavior across the package.
}

\item{level}{
Confidence level of the approximate confidence surfaces. The default is 0.95.
}
  \item{n.mix}{
 The number of simulations used to get the approximate confidence intervals or surfaces. If n.mix = 0, no simulation will be done and the face of the final projection will be used to compute the covariance matrix of the constrained estimate. The default is n.mix = 100L.
}
  \item{test}{
A logical scalar. If test == TRUE, then the p-value for the test \eqn{H_0:\theta} is in \eqn{V} versus \eqn{H_1:\theta} is in \eqn{C} is returned. \eqn{C} is the constraint cone of the form \eqn{\{\beta: A\beta \ge 0\}}, and \eqn{V} is the null space of \eqn{A}. The default is test = FALSE.
}


\item{subset}{
  Expression to select a subpopulation.
}

\item{...}{
  Extra arguments.
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
  If the response is Gaussian, type = "response" and type = "link" give the predicted mean; if the response is binomial, poisson or Gamma, type = "response" gives the predicted mean, and type = "link" gives the predicted systematic component.
}

\item{se.fit}{
Logical switch indicating if confidence intervals are required.
}

 
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
\item{amat}{The \eqn{k \times M} constraint matrix imposing shape constraints in each dimension, where \eqn{M} is the
total number of domains.}
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
\item{family}{The family parameter defined in a csvy formula.}
\item{df.residual}{The observed degree of freedom for the residuals of a csvy fit.}
\item{df.null}{The degree of freedom for the null model of a csvy fit.}
\item{domain}{Index of each domain in the data set contained in the survey.design object.}
\item{null.deviance}{The deviance for the null model of a csvy fit.}
\item{deviance}{The residual deviance of a csvy fit.}
  
%\item{terms}{TBA}
%\item{small_cells}{TBA}
%\item{prior.weights}{TBA}
%\item{grid_ps}{TBA}
%\item{obs_cells}{TBA}
%\item{edf}{TBA}
%\item{family}{TBA}
%\item{df.residual}{TBA}
%\item{df.null}{TBA}

%\item{ynm}{TBA}
%\item{fitted.values.un}{TBA}
%\item{domain}{TBA}
%\item{amat_0}{TBA}
%\item{data}{TBA}
%\item{weights}{TBA}
%\item{linear.predictors}{TBA}
%\item{fitted.values}{TBA}
%\item{xnms_add}{TBA}

%\item{w}{TBA}
%\item{y}{TBA}
%\item{bval}{TBA}
%\item{call}{TBA}
%\item{xmat_add0}{TBA}
%\item{empty_cells}{TBA}
%\item{zeros}{TBA}
%\item{cov.unscaled}{TBA}
%\item{ne}{TBA}

%\item{null.deviance}{TBA}
%\item{n.mix}{TBA}
%\item{linear.predictors.un}{TBA}
%\item{na.action}{TBA}

}

\details{
For binomial and Poisson families use \code{family=quasibinomial()}
and \code{family=quasipoisson()} to avoid a warning about non-integer
numbers of successes.  The `quasi' versions of the family objects give
the same point estimates and standard errors and do not give the
warning.

\code{predict} gives fitted values and sampling variability for specific new
values of covariates.  When \code{newdata} are the population mean it
gives the regression estimator of the mean, and when \code{newdata} are
the population totals and \code{total} is specified it gives the
regression estimator of the population total.  Regression estimators of
mean and total can also be obtained with \code{\link[survey]{calibrate}}.

%\code{coef}

%\code{vcov}

%\code{confint}

%\code{plotpersp}

%\code{summary}


}

\seealso{ 
  \code{\link[cgam]{plotpersp}}, to create a 3D Plot for a csvy Object
  
  \code{\link[cgam]{incr}}, to specify an increasing shape-restriction in a csvy Formula
  
  \code{\link[cgam]{decr}}, to specify an decreasing shape-restriction in a csvy Formula
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
# Example 1: monotonic in 1 predictor
data(nhdat2, package = 'csurvey')
#specify the design:
dstrat <- svydesign(ids = ~id, strata = ~str, data = nhdat2, weight = ~wt)

#uncomment to use parallel computing:
#options(csurvey.multicore=TRUE)
#mixture-variance-covariance matrix is simulated
set.seed(1)
ans <- csvy(chol ~ incr(age), design = dstrat, n.mix=5)

#check the constrained fit vs the unconstrained fit
summary(ans)
plot(ans, type = 'both')

\dontrun{
# Example 2: monotonic in 2 predictors and unconstrained in a 3rd predictor
data(nhdat2, package = 'csurvey')
#specify the design:
dstrat <- svydesign(ids = ~id, strata = ~str, data = nhdat2, weight = ~wt)

#use parallel computing:
#options(csurvey.multicore=TRUE)
#mixture-variance-covariance matrix is simulated
#Average cholestorel level increases in age, waist, and unconstrained in income
set.seed(1)
ans <- csvy(chol ~ incr(age)*incr(wcat)*icat, design = dstrat, test=FALSE, n.mix=5)

#visualize the constrained estimation with confidence bands
plot(ans, x1='icat', x2='wcat')

#create control object
ctl <- plot_csvy_control(ribbon_fill = "pink",x1lab = 'income',x2lab = 'waist')

plot(ans, x1='icat', x2='wcat', control=ctl)
}

# Example 3: example with a binomial response
library(NHANES)
library(survey)
data(NHANES)

nh <- subset(NHANES, !is.na(Education) & !is.na(BMI) & !is.na(Weight))
nh$DiabetesBin <- as.integer(nh$Diabetes == "Yes")
nh$BMIgroup <- cut(nh$BMI, breaks = c(0, 18.5, 25, 30, 35, 40, Inf), labels = FALSE)

# specify the design
dsgn <- svydesign(ids = ~1, strata = ~BMIgroup, weights = ~Weight, data = nh)

ans <- csvy(DiabetesBin ~ decr(Education) * incr(BMIgroup), design = dsgn,
  family = quasibinomial(link='logit'), n.mix=5)

summary(ans)

plot(ans, x1 = 'BMIgroup', x2 = 'Education')

ctl <- plot_csvy_control(
  x1size = 1.5, x2size = 2,
  angle = 45, hjust = .3)

plot(ans, x1 = 'Education', x2 = 'BMIgroup', control = ctl)


}

\keyword{main routine}
