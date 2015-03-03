\name{summary.clme}
\alias{summary.clme}

\title{ S3 method to summarize results for objects of class \code{clme} }
\description{ Summarizes the output of objects of class \code{clme}, such as those produced by \code{\link{clme}}. Prints a tabulated display of global and individual tests, as well as parameter estimates.}


\usage{
\method{summary}{clme}(object, alpha=0.05, digits=4, ...)
}

\arguments{
  \item{object}{an object of class \code{clme}. }
  \item{alpha}{level of significance.}
  \item{digits}{number of decimal digits to print.}
  \item{\dots}{ additional arguments passed to other functions. }
}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }


\note{ The individual tests are performed on the specified order. If no specific order was specified, then the individual tests are performed on the estimated order. }



\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}


\examples{

\dontrun{
data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 10)
                 
summary( clme.out )
}
}