\name{rat.blood}

\docType{data}

\alias{rat.blood}

\title{ Experiment on mice}

\description{
  This data set contains the data from an experiment on 24 Sprague-Dawley rats from Cora et al (2012).
}

\usage{rat.blood}
\format{A frame containing 241 observations on 13 variables.

 \tabular{rll}{
    [,1]  \tab id   \tab ID for rat (factor). \cr
    [,2]  \tab time \tab time period (in order, 0 , 6, 24, 48, 72, 96 hours). \cr
    [,3]  \tab temp \tab storage temperature reference (\code{''Ref''}) vs. room temperature (\code{''RT''}). \cr
    [,4]  \tab sex  \tab sex, male (\code{''Male''}) vs. female (\code{''Female''}). Coded as \code{''Female''=1}. \cr
    [,5]  \tab wbc  \tab white blood cell count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
    [,6]  \tab rbc  \tab red blood cell count )\eqn{10^6 / \mu L}{10^6 / mu L}). \cr
    [,7]  \tab hgb  \tab hemoglobin concentration (g/dl). \cr
    [,8]  \tab hct  \tab hematocrit (\%). \cr
    [,9]  \tab spun \tab (HCT \%). \cr
    [,10] \tab mcv  \tab MCV, a measurement of erythrocyte volume (fl). \cr
    [,11] \tab mch  \tab mean corpuscular hemoglobin (pg). \cr     % ????
    [,12] \tab mchc \tab mean corpuscular hemoglobin concentration (g/dl). \cr
    [,13] \tab plts \tab platelet count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
  }
}

\references{
Cora M, King D, Betz L, Wilson R, and Travlos G (2012). Artifactual changes in Sprauge-Dawley rat hematologic parameters after storage of samples at 3 C and 21 C. Journal of the American Association for Laboratory Animal Science, 51(5), 616-621. URL \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3447451/}.
 
}

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\keyword{datasets}