\name{model.matrix.clme}

\alias{model.matrix.clme}


\title{ Extracts the fixed-effects design matrix }
\description{
Extracts the fixed-effects design matrix from objects of class \code{clme}. 
}

\usage{
\method{model.matrix}{clme}(object, ... )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}


\value{
Returns a matrix.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

\dontrun{
data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 0)

model.matrix.clme( clme.out )
}

}

