\name{plotpersp_csvy_control}
\alias{plotpersp_csvy_control}
\title{Create a Control List for \code{plotpersp.csvy}}
\description{
Constructs a list of control parameters for use with \code{\link{plotpersp.csvy}}. 
It extends the default settings from \code{\link[cgam]{plotpersp_control}} with options specific to \code{csvy} plots.
}
\usage{
plotpersp_csvy_control(surface = c("C", "U"), x1nm = NULL, x2nm = NULL, categ = NULL,
  col = NULL, random = FALSE, ngrid = 12, xlim = NULL, ylim = NULL, zlim = NULL,
  xlab = NULL, ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL, main = NULL,
  sub = NULL, ticktype = "simple", ci = c("none", "lwr", "up", "both"), 
  cex = 1, categnm = NULL, type = c("response", "link"),
  cex.main = 0.8, box = TRUE, axes = TRUE, nticks = 5, palette = NULL, 
  NCOL = NULL, transpose = FALSE)
}
\arguments{
  \item{surface}{Plot the constrained (\code{"C"}) or unconstrained (\code{"U"}) fit.}
  \item{x1nm, x2nm}{Character strings naming the predictor variables for x and y axes.}
  \item{categ}{Optional character string naming a categorical covariate to stratify plots.}
  \item{col}{Color(s) for surfaces. Can be a palette or custom colors.}
  \item{random}{If \code{TRUE}, colors are assigned randomly.}
  \item{ngrid}{Number of grid points along each axis.}
  \item{xlim, ylim, zlim}{Optional limits for x, y, z axes.}
  \item{xlab, ylab, zlab}{Axis labels.}
  \item{th, ltheta}{Viewing and lighting angles for the plot.}
  \item{main, sub}{Plot title and subtitle.}
  \item{ticktype}{Type of ticks: \code{"simple"} or \code{"detailed"}.}
  \item{ci}{Confidence interval display mode: \code{"none"}, \code{"lwr"}, \code{"up"}, \code{"both"}.}
  \item{cex}{Scaling factor for labels.}
  \item{categnm}{Labels for each level of the categorical covariate.}
  \item{type}{Scale of the surface: \code{"response"} (default) or \code{"link"}.}
  \item{cex.main}{Scaling factor for main title text.}
  \item{box, axes}{Logical flags to show box and axes.}
  \item{nticks}{Number of tick marks along axes.}
  \item{palette}{Vector of colors for multi-surface plots.}
  \item{NCOL}{Number of columns in multi-panel layout.}
  \item{transpose}{Logical; if \code{TRUE}, transpose layout of multi-panel plots.}
}
\value{
A named list of graphical settings for use in \code{\link{plotpersp.csvy}}.
}
\seealso{
\code{\link{plotpersp.csvy}}, \code{\link[cgam]{plotpersp_control}}, \code{\link[graphics]{persp}}
}
\examples{
ctrl <- plotpersp_csvy_control(col = "topo", ci = "both", transpose = TRUE)
}
\keyword{graphics}
