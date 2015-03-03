\name{ranef.clme}

\alias{ranef.clme}
\alias{random.effects.clme}


\title{ Extract random effects }
\description{
Extracts the random effects estimates from objects of class \code{clme}. 
}

\usage{
\method{ranef}{clme}(object, ... )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}


\value{
Returns a numeric vector.
}


% \author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )
clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 0)

ranef.clme( clme.out )


}

