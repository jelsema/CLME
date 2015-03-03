\name{AIC.clme}

\alias{AIC.clme}


\title{Akaike information criterion }
\description{
Calculates the Akaike and Bayesian information criterion for objects of class \code{clme}. 
}

\usage{
\method{AIC}{clme}(object, ..., k=2 )
}

\arguments{
  \item{object}{ object of class \code{\link{clme}}. }
  \item{...}{space for additional arguments.}  
  \item{k}{ value multiplied by number of coefficients}
}


\details{
The log-likelihood is assumed to be the Normal distribution. The model uses residual bootstrap methodology, and Normality is neither required nor assumed. Therefore the log-likelihood and these information criterion may not be useful measures for comparing models.

For \code{k=2}, the function computes the AIC. To obtain BIC, set \eqn{k = log( n/(2*pi) )}; which the method \code{BIC.clme} does.
}

\value{
Returns the information criterion (numeric).
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

AIC( clme.out )
AIC( clme.out, k=log( nobs(clme.out)/(2*pi) ) )

}

