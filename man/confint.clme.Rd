\name{confint.clme}

\alias{confint.clme}


\title{ Individual confidence intervals }
\description{
Calculates confidence intervals for fixed effects parameter estimates in objects of class \code{clme}. 
}

\usage{
\method{confint}{clme}(object, parm, level=0.95, ... )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{parm}{parameter for which confidence intervals are computed (not used).}
  \item{level}{nominal confidence level.}
  \item{...}{space for additional arguments.}  
}


\details{
Confidence intervals are computed using Standard Normal critical values. Standard errors are taken from the covariance matrix of the unconstrained parameter estimates.
}

\value{
Returns a matrix with two columns named lcl and ucl (lower and upper confidence limit).
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

confint.clme( clme.out )


}

