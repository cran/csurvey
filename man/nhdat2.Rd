\name{nhdat2}
\alias{nhdat2}

\title{
  A Subset of National Health and Nutrition Examination Survey (NHANES) 
}

\description{

  The National Health and Nutrition Examination Survey (NHANES) combines in-person interviews and physical examinations to produce a comprehensive data set from a probability sample of residents of the U.S.  
  
  This data set is a subset of the NHANES data with 1,933 subjects.

}

\usage{data(nhdat2)}
\format{
  A data frame with 1933 observations on the following 8 variables.
  \describe{
   \item{\code{id}}{a identification vector specifying cluster ids from largest level to smallest level}
    \item{\code{chol}}{a continuous vector of cholesterol level}
    \item{\code{wcat}}{a vector of categorized waist and height ratio}
    \item{\code{icat}}{an ordinal vector of categorized income level}
    \item{\code{gender}}{a binary vector of genders}
    \item{\code{age}}{a vector of categorized age}
    \item{\code{wt}}{sampling weight within each stratum}
    \item{\code{str}}{a numeric vector}
  }
}
\details{
  The variable chol in this data set is continuous, which the varible chol in the nhdat data set is binomial.
}

\examples{
\dontrun{
data(nhdat2)
summary(nhdat2)
}
}

\keyword{datasets}
