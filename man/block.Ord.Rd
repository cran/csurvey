\name{block.Ord}

\alias{block.Ord}

\title{Specify a Block Monotonic Shape-Restriction in a CSVY Formula}

\description{
  A symbolic routine to define that a vector of domain means follows a monotonic ordering in a predictor in a formula argument to csvy. This is the unsmoothed version.
}

\usage{
  block.Ord(x, order = NULL, numknots = 0, knots = 0, space = "E")
}

\arguments{
 \item{x}{A numeric predictor which has the same length as the response vector.}

  \item{order}{
A \eqn{1 \times M} vector defining the order of domains when the shape constraint is block ordering.
}
   \item{numknots}{The number of knots used to smoothly constrain a predictor. The value should be \eqn{0} for a shape-restricted predictor without smoothing. The default value is \eqn{0}.}

  \item{knots}{The knots used to smoothly constrain a predictor. The value should be \eqn{0} for a shape-restricted predictor without smoothing. The default value is \eqn{0}.}

 \item{space}{A character specifying the method to create knots. It will not be used for a shape-restricted predictor without smoothing. The default value is "E".}

}

\value{The vector x with five attributes, i.e., name: the name of x; shape: 9("block ordering"); numknots: the numknots argument in "block.Ord"; knots: the knots argument in "block.Ord"; space: the space argument in "block.Ord".}


\author{
  Xiyue Liao
}

\seealso{
  \code{\link{incr}}, \code{\link{decr}}
}

\keyword{shape routine}

