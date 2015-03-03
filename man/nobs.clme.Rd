\name{nobs.clme}

\alias{nobs.clme}


\title{ Number of observations }
\description{
Obtains the number of observations used to fit an model for objects of class \code{clme}. 
}

\usage{
\method{nobs}{clme}(object, ... )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}


\value{
Numeric.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 0)

nobs.clme( clme.out )

}

