\name{clme_resids}

\alias{clme_resids}


\title{ Computes various types of residuals }
\description{
Computes several types of residuals for objects of class \code{clme}. 
}

\usage{
clme_resids( formula, data, gfix=NULL, ncon=1 )
}

\arguments{
  \item{formula}{ a formula expression. The constrained effect(s) must come before any unconstrained covariates on the right-hand side of the expression. The first \code{ncon} terms will be assumed to be constrained.}
  \item{data}{ data frame containing the variables in the model. }
  \item{gfix}{ optional vector of group levels for residual variances. Data should be sorted by this value.}
  \item{ncon}{ the number of variables in \code{formula} that are constrained.}
}

\details{
For fixed-effects models \eqn{Y = X\beta + \epsilon}{Y = X*b + e}, residuals are given as \eqn{\hat{e} = Y - X\hat{\beta}}{ ehat = Y - X*betahat}.

For mixed-effects models \eqn{Y = X\beta + + U\xi + \epsilon}{Y = X*b + U*xi + e}, three types of residuals are available.
\eqn{PA = Y - X\hat{\beta}}{ PA = Y - X*betahat}\\
\eqn{SS = U\hat{\xi}}{ SS = U*xihat}\\
\eqn{FM = Y - X\hat{\beta} - U\hat{\xi}}{ FM = Y - X*betahat - U*xihat}
}

\value{
List containing the elements \code{PA}, \code{SS}, \code{FM}, \code{cov.theta}, \code{xi}, \code{ssq}, \code{tsq}.

\code{PA}, \code{SS}, \code{FM} are defined above (for fixed-effects models, the residuals are only \code{PA}). Then \code{cov.theta} is the unconstrained covariance matrix of the fixed-effects coefficients, \code{xi} is the vector of random effect estimates, and \code{ssq} and \code{tsq} are unconstrained estimates of the variance components.

}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

\dontrun{

data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )

clme.out <- clme_resids(mcv ~ time + temp + sex + (1|id), data = rat.blood, ncon=1 )

}

}

