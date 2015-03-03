\name{is.clme}

\alias{as.clme}
\alias{is.clme}


\title{Constructor method for objects S3 class clme}
\description{
Test if an object is of class \code{clme} or coerce an object to be such.
}

\usage{
\method{as}{clme}(x, ... )
\method{is}{clme}( x )
}

\arguments{
  \item{x}{ list with the elements corresponding to the output of \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}


\value{
Returns an object of the class \code{clme}.
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
                 
is.clme( clme.out )
as.clme( clme.out )

}

