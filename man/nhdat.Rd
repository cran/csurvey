\name{nhdat}
\docType{data}
\alias{nhdat}
\title{Subset of NHANES 2009 to 2010 Data for Constrained Survey Analysis}

\description{
A structured subset of the 2009 to 2010 NHANES data designed for illustrating constrained survey estimation methods.
}

\usage{
data(nhdat)
}

\format{
A data frame with 1,680 observations on the following 7 variables:
\describe{
  \item{\code{id}}{Cluster identifier derived from NHANES sequence number (\code{SEQN}).}
  \item{\code{chol}}{Binary indicator of high total cholesterol: 1 if total cholesterol > 200 mg/dL, 0 otherwise. Derived from \code{LBXTC} in \code{TCHOL_F.XPT}.}
  \item{\code{wcat}}{Four-level ordinal variable for waist to height ratio category, based on \code{BMXWAIST} and \code{BMXHT} from \code{BMX_F.XPT}.}
  \item{\code{gender}}{Gender of the participant (1 = male, 2 = female), from \code{RIAGENDR} in \code{DEMO_F.XPT}.}
  \item{\code{age}}{Age in years (continuous), from \code{RIDAGEYR} in \code{DEMO_F.XPT}.}
  \item{\code{wt}}{Sampling weight within strata, based on \code{WTINT2YR} from \code{DEMO_G.XPT}.}
  \item{\code{str}}{Stratum identifier, based on \code{SDMVPSU} from \code{DEMO_G.XPT}.}
}
}

\details{
This subset includes participants aged 21 to 45 years, selected for illustrating the estimation of the probability of high cholesterol using order constrained survey methods.
}

\source{
National Center for Health Statistics. NHANES 2009 to 2010 Public Use Data Files.\\
\url{https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2009}
}

%\keywords{datasets}
