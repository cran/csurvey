\name{incr}

\alias{incr}

\title{Specify an Increasing Shape-Restriction in a CSVY Formula}

\description{
A symbolic routine to define that a vector of domain means is increasing in a predictor in a formula argument to csvy. This is the unsmoothed version.}

\usage{incr(x, numknots = 0, knots = 0, space = "E")}

\arguments{
  \item{x}{A numeric predictor which has the same length as the response vector.}
  \item{numknots}{The number of knots used to smoothly constrain a predictor. The value should be \eqn{0} for a shape-restricted predictor without smoothing. The default value is \eqn{0}.}
  \item{knots}{The knots used to smoothly constrain a predictor. The value should be \eqn{0} for a shape-restricted predictor without smoothing. The default value is \eqn{0}.}
 \item{space}{A character specifying the method to create knots. It will not be used for a shape-restricted predictor without smoothing. The default value is "E".}
}

\value{The vector x with five attributes, i.e., name: the name of x; shape: 1("increasing"); numknots: the numknots argument in "incr"; knots: the knots argument in "incr"; space: the space argument in "incr".}

\details{
 The subroutine "incr" is borrowed from the \emph{R} package \emph{cgam}. For now, the arguments numknots, knots, and space will not be used because only unsmoothed shape-constrained estimation is provided in this package.

  "incr" is a helper function to create the constraint matrix used in csvy.

  See references cited in this section for more details.
}

\author{
  Xiyue Liao
}

\references{
  Liao, X. and Meyer, M. C. (2019) cgam: An R Package for the Constrained Generalized Additive Model. \emph{Journal of Statistical Software \bold{89(5)}}.
}

\seealso{
  \code{\link[cgam]{incr}}, \code{\link{csvy}}
}

\keyword{shape routine}
