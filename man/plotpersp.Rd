\name{plotpersp}

\alias{plotpersp}

\title{Create a 3D Plot for a CSVY Object}

\description{This routine is borrowed from the cgam package. Given an object of the csvy class, which has at least two non-parametrically modelled predictors, this routine will make a 3D plot of the estimated domain mean surface with a set of two non-parametrically modelled predictors in the formula being the x and y labs. 

This routine is an extension of the generic R graphics routine persp. 
}

\usage{
plotpersp(object, x1 = NULL, x2 = NULL,...)
}

\arguments{
  \item{object}{An object of the csvy class with at least two non-parametrically modelled predictors.}
    \item{x1}{A non-parametrically modelled predictor in a csvy fit. If the user omits x1 and x2, then the first two non-parametric predictors in a csvy formula will be used.}
  \item{x2}{A non-parametrically modelled predictor in a csvy fit. If the user omits x1 and x2, then the first two non-parametric predictors in a csvy formula will be used.}
\item{...}{
Arguments to be passed to the S3 method for the csvy class:
\itemize{
  \item{x1nm: }{Character name of x1.}
  \item{x2nm: }{Character name of x2.}
  %\item{data: }{The data frame based on which the user get a csvy fit.}
  \item{col: }{The color(s) of a 3D plot created by plotpersp. If col == NULL, "white" will be used when there is only one surface in the plot. The default is col =  NULL.}
  \item{xlim: }{The xlim argument inherited from the persp routine.}
  \item{ylim: }{The ylim argument inherited from the persp routine.}
  \item{zlim: }{The zlim argument inherited from the persp routine.}
  \item{xlab: }{The xlab argument inherited from the persp routine.}
  \item{ylab: }{The ylab argument inherited from the persp routine.}
  \item{zlab: }{The zlab argument inherited from the persp routine.}
  \item{main: }{The main argument inherited from the persp routine.}
 \item{categ: }{If there are more than 2 predictors in the csvy fit, then the user can specify the argument categ to be a character representing a third predictor \eqn{x_3} in the formula, then multiple 3D plots, which represent the domains of \eqn{x_3} in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. The default is categ = NULL.}
  \item{th: }{The theta argument inherited from the persp routine.}
  \item{ltheta: }{The ltheta argument inherited from the persp routine.}
  \item{main: }{The main argument inherited from the persp routine.}
  \item{ticktype: }{The ticktype argument inherited from the persp routine.}
  %\item{\dots}{Further arguments passed to the routine.}
}
}
}


\details{
  The graphic routine "plotpersp" is borrowed from the \emph{R} package \emph{cgam}.

  See the help page of \code{\link{csvy}} for examples.
}

\value{The routine plotpersp returns a 3D plot of an object of the csvy class. The \eqn{x} lab and \eqn{y} lab represent a set of non-parametrically modelled predictors used in a csvy formula, and the \eqn{z} lab represents estimated domain means.}


\references{The official documentation for the generic R routine persp: \url{http://stat.ethz.ch/R-manual/R-patched/library/graphics/html/persp.html}}

\author{
  Xiyue Liao
}

\seealso{
   \code{\link[cgam]{plotpersp}}, \code{\link{csvy}}
}


\keyword{3D plot routine}
