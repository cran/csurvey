\name{nhdat}
\alias{nhdat}
\title{
  A Subset of National Health and Nutrition Examination Survey (NHANES) 
}

\description{

  The National Health and Nutrition Examination Survey (NHANES) combines in-person interviews and physical examinations to produce a comprehensive data set from a probability sample of residents of the U.S.  
  
  This data set is a subset of the NHANES data with 1,680 subjects.

}
\usage{data(nhdat)}
\format{
  \describe{
    \item{\code{id}}{a identification vector specifying cluster ids from largest level to smallest level}
    \item{\code{chol}}{a binomial vector showing cholestoral level. 1: high; 0: low}
    \item{\code{wcat}}{a vector of categorized waist and height ratio}
    \item{\code{gender}}{a binary vector of genders}
    \item{\code{age}}{a vector of categorized age}
    \item{\code{wt}}{sampling weight within each stratum}
    \item{\code{str}}{a numeric vector}
  }
}

\examples{
\dontrun{
data(nhdat)
summary(nhdat)
}
}

\keyword{datasets}
