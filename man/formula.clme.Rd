\name{formula.clme}

\alias{formula.clme}


\title{ Extract formula }
\description{
Extracts the formula from objects of class \code{clme}. 
}

\usage{
\method{formula}{clme}(x, ... )
}

\arguments{
  \item{x}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}

\details{ The package \pkg{CLME} parametrizes the model with no intercept term. If an intercept was included, it will be removed automatically. }

\value{
Returns a formula object.
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

formula.clme( clme.out )


}

