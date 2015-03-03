\name{VarCorr.clme}

\alias{VarCorr.clme}


\title{ Variance components }
\description{
Extracts variance components for objects of class \code{clme}. 
}

\usage{
\method{VarCorr}{clme}(object, ... )
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

VarCorr.clme( clme.out )


}

