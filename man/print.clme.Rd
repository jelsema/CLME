\name{print.clme}

\alias{print.clme}


\title{ Printout for fitted object }
\description{
Prints basic information on a fitted object of class \code{clme}. 
}

\usage{
\method{print}{clme}(x, ... )
}

\arguments{
  \item{x}{ object of class \code{\link{clme}} to be printed. }
  \item{...}{space for additional arguments.}  
}


\value{
Text printed to console.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

\dontrun{
set.seed( 42 )
data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 10)

print.clme( clme.out )
}

}

