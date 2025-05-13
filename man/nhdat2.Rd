\name{nhdat2}
\docType{data}
\alias{nhdat2}
\title{Subset of NHANES 2009 to 2010 Data for Constrained Survey Analysis}
\description{
A structured subset of the 2009 to 2010 NHANES data designed for illustrating constrained survey estimation methods.
}
\usage{data(nhdat2)}
\format{
A data frame with 1,933 observations and the following 8 variables:
\describe{
  \item{\code{id}}{An identification vector specifying cluster ids from largest level to smallest level, derived from NHANES sequence number \code{SEQN}.}
  \item{\code{chol}}{Total cholesterol, measured in mg/dL. This variable is derived from \code{LBXTC} in the laboratory file \code{TCHOL_F.XPT}.}
  \item{\code{wcat}}{A 4 level ordinal categorical variable representing waist to height ratio categories, computed from \code{BMXWAIST} (waist circumference in cm) and \code{BMXHT} (height in cm) in the body measures file \code{BMX_F.XPT}.}
  \item{\code{icat}}{A 4 level ordinal categorical variable. It categorizes income based on the ratio of family income to the federal poverty level (INDFMPIR), 
    with cutpoints at 0.75, 1.38, and 3.5 to reflect meaningful policy thresholds. It is derived from \code{INDFMIN2} in the demographics file \code{DEMO_F.XPT}.}
  \item{\code{gender}}{Gender of the participant, with values 1 (male) and 2 (female), derived from \code{RIAGENDR} in \code{DEMO_F.XPT}.}
  \item{\code{age}}{Age in years (continuous), derived from \code{RIDAGEYR} in the demographics file \code{DEMO_F.XPT}.}
  \item{\code{wt}}{Sampling weight within each stratum, derived from (\code{WTINT2YR}) from \code{DEMO_G.XPT}.}
  \item{\code{str}}{Stratum identifier, derived from (\code{SDMVPSU}) from \code{DEMO_G.XPT}.}
}
}
\details{
  This subset includes participants aged 21 through 45, selected to demonstrate estimation of domain means using order constrained methods. 
}
\source{
  National Center for Health Statistics. NHANES 2009 to 2010 Public Use Data Files.  
  \url{https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2009}
}
%\keywords{datasets}