\name{fibroid}

\docType{data}

\alias{fibroid}

\title{ Fibroid Growth Study }

\description{
  This data set contains a subset of the data from the Fibroid Growth Study.
}

\usage{fibroid}
\format{A frame containing 240 observations on 9 variables.

 \tabular{rll}{
    [,1] \tab  ID  \tab ID for subject. \cr
    [,2] \tab fid  \tab ID for fibroid (each women could have multiple fibroids). \cr
    [,3] \tab lfgr \tab log fibroid growth rate. See details. \cr
    [,4] \tab age  \tab age category Younger, Middle, Older. \cr
    [,5] \tab loc  \tab location of fibroid, corpus, fundus, or lower segment. \cr
    [,6] \tab bmi  \tab body mass index of subject. \cr
    [,7] \tab preg \tab parity, whether the subject had delivered a child. \cr
    [,8] \tab race \tab race of subject (Black or White only). \cr
    [,9] \tab vol  \tab initial volume of fibroid. \cr
  }
}


\details{

The response variable \code{lfgr} was calculated as the change in log fibroid volume, divided by the length of time between measurements. The growth rates were averaged to produce a single value for each fibroid, which was scaled to represent a 6-month percent change in volume.

}

\references{
Peddada, Laughlin, Miner, Guyon, Haneke, Vahdat, Semelka, Kowalik, Armao, Davis, and Baird(2008). Growth of Uterine Leiomyomata Among Premenopausal Black and White Women. Proceedings of the National Academy of Sciences of the United States of America, 105(50), 19887-19892. URL \url{http://www.pnas.org/content/105/50/19887.full.pdf}.
}




\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\keyword{datasets}