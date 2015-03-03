\name{print.varcorr_clme}

\alias{print.varcorr_clme}


\title{ Printout for variance components }
\description{
Prints variance components of an objects of \code{clme}. 
}

\usage{
\method{print}{varcorr_clme}(x, ... )
}

\arguments{
  \item{x}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}


\value{
Text printed to console.
}


% \author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

\dontrun{
data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 0)

print.varcorr_clme( clme.out )
}

}

