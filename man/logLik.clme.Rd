\name{logLik.clme}

\alias{logLik.clme}


\title{ Log-likelihood }
\description{
Computes the log-likelihood of the fitted model for objects of class \code{clme}. 
}

\usage{
\method{logLik}{clme}(object, ... )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
}

\details{ The log-likelihood is computed using the Normal distribution. The model uses residual bootstrap methodology, and Normality is neither required nor assumed. Therefore the log-likelihood may not be useful measures in the context of \pkg{CLME}.
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

logLik.clme( clme.out )


}

